FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y build-essential

COPY environment.yaml /environment.yaml
RUN conda env create -f /environment.yaml

ENV PATH /opt/conda/envs/transcriptomics/bin:$PATH

WORKDIR /pipeline
COPY . /pipeline

CMD ["bash", "run_pipeline.sh"]
