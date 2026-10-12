#!/usr/bin/env python3
"""Compile in an isolated directory and render every page for visual review."""
from pathlib import Path
import hashlib
import io
import json
import os
import shutil
import subprocess
from PIL import Image, ImageDraw
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / 'reduced_search/deliverables/Grothendieck_reduced_sequence_2026-10-01.tex'
BUILD = ROOT / 'tonight/latex_final_build'
QA = ROOT / 'tonight/pdf_qa/final'

def atomic_bytes(path, data):
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('wb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)

def main():
    if BUILD.exists(): shutil.rmtree(BUILD)
    if QA.exists(): shutil.rmtree(QA)
    BUILD.mkdir(parents=True)
    QA.mkdir(parents=True)
    shutil.copyfile(SOURCE, BUILD / SOURCE.name)
    for iteration in (1, 2):
        run = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                              SOURCE.name], cwd=BUILD, capture_output=True, check=True)
        (BUILD / f'pass{iteration}.txt').write_bytes(run.stdout + run.stderr)
    log = (BUILD / SOURCE.with_suffix('.log').name).read_text()
    assert not any(word in log for word in ('LaTeX Warning:', 'Overfull', 'Underfull'))
    pdf = BUILD / SOURCE.with_suffix('.pdf').name
    reader = PdfReader(pdf)
    assert all(len(page.extract_text()) > 100 for page in reader.pages)
    atomic_bytes(SOURCE.with_suffix('.pdf'), pdf.read_bytes())
    thumbs = []
    for index in range(len(reader.pages)):
        run = subprocess.run(['pdftoppm', '-f', str(index+1), '-l', str(index+1),
                              '-scale-to', '1600', '-png', str(pdf)],
                             capture_output=True, check=True)
        page = Image.open(io.BytesIO(run.stdout)); page.load()
        assert page.width > 500 and page.height > 500
        atomic_bytes(QA / f'page-{index+1:02d}.png', run.stdout)
        thumb = page.convert('RGB'); thumb.thumbnail((320, 414))
        panel = Image.new('RGB', (330, 438), 'white')
        panel.paste(thumb, ((330-thumb.width)//2, 22))
        ImageDraw.Draw(panel).text((6, 4), f'Page {index+1}', fill='black')
        thumbs.append(panel)
    for start in range(0, len(thumbs), 9):
        chosen = thumbs[start:start+9]
        sheet = Image.new('RGB', (990, 438*((len(chosen)+2)//3)), '#dddddd')
        for offset, panel in enumerate(chosen): sheet.paste(panel, ((offset%3)*330, (offset//3)*438))
        stream = io.BytesIO(); sheet.save(stream, format='PNG')
        atomic_bytes(QA / f'contact-{start+1:02d}.png', stream.getvalue())
    result = dict(status='PASS_COMPILE_PARSE_AND_RENDER_ALL_PAGES_REQUIRES_VISUAL_REVIEW',
                  pages=len(reader.pages), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                  pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
                  pdf_bytes=pdf.stat().st_size, qa_directory=str(QA))
    (QA / 'BUILD_CHECK.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__': main()
