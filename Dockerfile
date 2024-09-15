# --- build layer
FROM python:3.12-bookworm AS builder

# TODO look into removing this dep
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libopenblas-dev \
    libxrender-dev \
    libxext-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install poetry==1.8.3
RUN poetry self add poetry-plugin-bundle

ENV VIRTUAL_ENV=/app/.venv \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1

WORKDIR /app

COPY src ./src
COPY README.md pyproject.toml poetry.lock ./
RUN touch README.md

RUN poetry bundle venv $VIRTUAL_ENV \
    --without dev \
    --compile \
    --verbose \
    --no-cache

# --- runtime layer
FROM python:3.12-slim-bookworm AS runtime

# FIXME pixi using conda bin deps should fix
# potentially other ways of fixing this as well
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libopenblas-dev \
    libxrender-dev \
    libxext-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ARG NB_USER=jovyan
ARG NB_UID=1000

ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH" \
    USER=${NB_USER} \
    NB_UID=${NB_UID} \
    HOME=/home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

COPY --from=builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}

# copy contents and change their ownership
COPY artifacts/ ${HOME}/artifacts
COPY notebooks/ ${HOME}/notebooks
COPY voila.json ${HOME}/voila.json

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR ${HOME}

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser"]
