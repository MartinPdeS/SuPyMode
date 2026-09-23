.DEFAULT_GOAL := help
PYTHON ?= python3
BUILD_DIR ?= build
ROOT_DIR := $(CURDIR)
TAG_VERSION ?= $(or $(VERSION),$(filter v%,$(MAKECMDGOALS)))
PYBIND11_DIR = $(shell $(PYTHON) -m pybind11 --cmakedir)
RELEASE_KIND := $(filter major minor patch,$(MAKECMDGOALS))

.PHONY: help setup configure build install editable test quality check release-check tag release major minor patch clean

ifneq ($(filter tag,$(MAKECMDGOALS)),)
ifneq ($(strip $(TAG_VERSION)),)
.PHONY: $(TAG_VERSION)
$(TAG_VERSION):
	@:
endif
endif

help:
	@echo "SuPyMode development commands"
	@echo "  make setup                 Install editable development dependencies"
	@echo "  make editable              Build and install an editable package"
	@echo "  make test                  Run the test suite"
	@echo "  make quality               Run static checks"
	@echo "  make release-check         Check version metadata consistency"
	@echo "  make tag VERSION=vX.Y.Z    Create a release commit and annotated tag"
	@echo "  make release patch         Create and push the next patch release"
	@echo "  make release minor         Create and push the next minor release"
	@echo "  make release major         Create and push the next major release"

setup:
	$(PYTHON) -m pip install -e ".[testing,documentation,dev]"
configure:
	cmake -S . -B $(BUILD_DIR) -Dpybind11_DIR="$(PYBIND11_DIR)" -DPython_EXECUTABLE="$$(which $(PYTHON))" -DCMAKE_INSTALL_PREFIX="$(ROOT_DIR)"
build:
	cmake --build $(BUILD_DIR) -j
install:
	cmake --install $(BUILD_DIR)
editable:
	$(PYTHON) -m pip install --no-build-isolation -e .
test:
	$(PYTHON) -m pytest --config-file=pytest.ini
quality:
	$(PYTHON) -m ruff check SuPyMode tests
check: quality test
release-check:
	$(PYTHON) tools/check_release.py $(if $(VERSION),--version $(VERSION),)
tag:
	$(PYTHON) tools/release_tag.py "$(TAG_VERSION)"
release:
	@test "$(words $(RELEASE_KIND))" -eq 1 || { echo "usage: make release [patch|minor|major]" >&2; exit 2; }
	@set -eu; release_tag="$$($(PYTHON) tools/next_release_version.py $(RELEASE_KIND))"; \
	$(PYTHON) tools/release_tag.py "$$release_tag"; \
	git push origin HEAD "refs/tags/$$release_tag"
major minor patch:
	@:
clean:
	rm -rf $(BUILD_DIR)
