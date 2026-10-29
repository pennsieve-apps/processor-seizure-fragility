.PHONY: help image test test-local fixture run clean

SERVICE_NAME ?= processor-seizure-fragility
IMAGE ?= $(SERVICE_NAME):dev
PYTHON ?= python3
export DOCKER_DEFAULT_PLATFORM = linux/amd64

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test        - unit tests (tests/testthat, synthetic EDF fixtures written in R) inside the image"
	@echo "make test-local  - the same tests with a local R that has EZFragility 2.1.1, gsignal, jsonlite, testthat"
	@echo "make fixture     - synthetic BIDS-iEEG derivative + seizure_onsets.csv -> data/fixture/{ieeg,onsets}"
	@echo "make run         - build the image and run it over data/fixture/, writing to data/output/"
	@echo "make clean       - remove data/output/ and data/fixture/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (make fixture needs NumPy only)"

image:
	docker build --platform linux/amd64 -t $(IMAGE) .

test: image
	docker run --rm --platform linux/amd64 --memory 2g --entrypoint Rscript $(IMAGE) /app/tests/testthat.R

test-local:
	Rscript tests/testthat.R

fixture:
	$(PYTHON) tests/make_fixture.py data/fixture

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort
	@head -5 data/output/fragility.csv

clean:
	rm -rf data/output data/fixture
