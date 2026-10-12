#!/usr/bin/env python3
"""Resolve original logical labels through a validated compact-package map."""
from pathlib import Path
import json
import re


def _unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise RuntimeError('Duplicate paths.json key: ' + str(key))
        answer[key] = value
    return answer


def canonical(name, *, allow_parents=False):
    """Normalize archived aliases, refusing absolute paths and root escapes."""
    if not isinstance(name, str) or not name or '\\' in name or ':' in name or '\x00' in name:
        raise RuntimeError('Invalid relative package path: ' + repr(name))
    if name.startswith('/'):
        raise RuntimeError('Absolute package path: ' + name)
    parts = []
    for part in name.split('/'):
        if part in ('', '.'):
            if not allow_parents:
                raise RuntimeError('Noncanonical package path: ' + name)
            continue
        if part == '..':
            if not allow_parents or not parts:
                raise RuntimeError('Package path traversal: ' + name)
            parts.pop()
        else:
            parts.append(part)
    if not parts:
        raise RuntimeError('Empty package path: ' + name)
    return '/'.join(parts)


class PackagePaths:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.map_path = self.root / 'paths.json'
        data = json.loads(self.map_path.read_text(), object_pairs_hook=_unique_object)
        if data.get('format') != 1 or not isinstance(data.get('files'), dict):
            raise RuntimeError('Unsupported paths.json format')
        self.files = data['files']
        physical = {}
        prefixes = {}
        logical = {}
        for name, target in self.files.items():
            if canonical(name) != name or canonical(target) != target:
                raise RuntimeError('Noncanonical paths.json member: ' + name)
            folded_name = name.casefold()
            if folded_name in logical:
                raise RuntimeError('Logical case collision: ' + name)
            logical[folded_name] = name
            folded = target.casefold()
            if folded in physical:
                raise RuntimeError('Physical target or case collision: ' + target)
            physical[folded] = target
            parts = target.split('/')
            for part in parts:
                if (len(part) > 255 or part.endswith(('.', ' '))
                        or any(ord(ch) < 32 or ch in '<>"|?*' for ch in part)
                        or re.fullmatch(r'(?:CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\..*)?', part, re.I)):
                    raise RuntimeError('Windows-incompatible physical member: ' + target)
            for count in range(1, len(parts)):
                prefix = '/'.join(parts[:count])
                old = prefixes.setdefault(prefix.casefold(), prefix)
                if old != prefix:
                    raise RuntimeError('Physical directory case collision: ' + target)
        if physical.keys() & prefixes.keys():
            raise RuntimeError('Physical file/directory collision')
        for target in self.files.values():
            path = self.root / target
            if not path.is_file() or not path.resolve().is_relative_to(self.root):
                raise RuntimeError('Missing or escaping mapped member: ' + target)

    def resolve(self, logical, *, prefix=None):
        name = canonical(logical, allow_parents=True)
        if prefix is not None:
            prefix = canonical(prefix)
            if not name.startswith(prefix + '/'):
                raise RuntimeError('Logical path leaves required root: ' + logical)
        try:
            return self.root / self.files[name]
        except KeyError:
            raise RuntimeError('Unmapped logical package member: ' + name) from None

    def frozen(self, name):
        return self.resolve('verification/frozen/' + name, prefix='verification/frozen')
