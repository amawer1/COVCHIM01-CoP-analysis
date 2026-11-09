.PHONY: install_mac install_linux install_windows install_r_packages

# macOS: system libraries via Homebrew, then R packages
install_mac:
	brew install nlopt
	$(MAKE) install_r_packages

# Debian/Ubuntu: system libraries needed to build the R packages from source, then R packages
install_linux:
	sudo apt-get update
	sudo apt-get install -y \
		libnlopt-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
		libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
		libpng-dev libtiff5-dev libjpeg-dev libcairo2-dev
	$(MAKE) install_r_packages

# Windows: CRAN provides pre-built binaries, so no system libraries are needed
install_windows: install_r_packages

install_r_packages:
	Rscript install_packages.R
