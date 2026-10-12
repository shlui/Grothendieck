FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends libmpfr6 libgmp10 && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONINTMAXSTRDIGITS=0
WORKDIR /repo
COPY . /repo
ENTRYPOINT ["python", "tools/replay.py"]
