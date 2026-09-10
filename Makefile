# Docker-based development environment for mapfit.
# The host does not need R installed; everything runs inside the image
# built from ./Dockerfile. See that file for the available targets.

IMAGE := mapfit-dev

# Run as the invoking user so that generated files (man/, NAMESPACE,
# mapfit.Rcheck/) are not left owned by root on the host.
# HOME must point somewhere writable, otherwise R fails to start.
DOCKER_RUN := docker run --rm \
	-v "$(CURDIR)":/pkg \
	-w /pkg \
	-u $(shell id -u):$(shell id -g) \
	-e HOME=/tmp \
	$(IMAGE)

.PHONY: image document test check win mac submit readme shell clean

image:
	docker build -t $(IMAGE) .

document:
	$(DOCKER_RUN) Rscript -e 'roxygen2::roxygenise()'

test:
	$(DOCKER_RUN) Rscript -e 'devtools::test()'

check:
	$(DOCKER_RUN) sh -c 'R CMD build . && R CMD check --as-cran mapfit_*.tar.gz'

# Pre-submission checks on the platforms this image cannot provide.
# Both upload the tarball to a build service and mail the result to the
# Maintainer address in DESCRIPTION -- nothing useful appears on stdout.
win:
	$(DOCKER_RUN) Rscript -e 'devtools::check_win_devel()'

mac:
	$(DOCKER_RUN) Rscript -e 'devtools::check_mac_release()'

# Submit to CRAN. Uploads the tarball together with cran-comments.md as the
# "Optional comment", then records CRAN-SUBMISSION. CRAN still mails a
# confirmation link that a human has to click.
#
# Open an interactive R session in the package, from which
# devtools::submit_cran() uploads the tarball with cran-comments.md as the
# submission comment and writes CRAN-SUBMISSION.
#
# The call is not run automatically: submit_cran() confirms through
# utils::menu(), which reads the terminal, so it cannot be passed with -e
# ("R -e" feeds the expression through stdin, leaving menu() nothing to read;
# "R --interactive -e" ignores the expression). Running it from .Rprofile or
# .First does not work either -- utils is not attached that early.
submit:
	@echo ''
	@echo '  In the R session that opens, type:'
	@echo ''
	@echo '      devtools::submit_cran()'
	@echo ''
	@echo '  Answer the two confirmations, then click the link CRAN mails you.'
	@echo '  q() leaves without submitting.'
	@echo ''
	docker run --rm -it \
		-v "$(CURDIR)":/pkg \
		-w /pkg \
		-u $(shell id -u):$(shell id -g) \
		-e HOME=/tmp \
		$(IMAGE) R --quiet --no-save

readme:
	$(DOCKER_RUN) Rscript -e 'devtools::build_readme()'

shell:
	docker run --rm -it \
		-v "$(CURDIR)":/pkg \
		-w /pkg \
		-u $(shell id -u):$(shell id -g) \
		-e HOME=/tmp \
		$(IMAGE) bash

clean:
	rm -rf mapfit.Rcheck mapfit_*.tar.gz src/*.o src/*.so
