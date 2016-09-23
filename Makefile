KB_TOP ?= /kb/dev_container
TARGET ?= /kb/deployment
MODULE = RNASeq
MODULE_CAPS = KBaseRNASeq
#MODULE_PORT = 6606
SPEC_FILE = KBaseRNASeq.spec

#End of user defined variables

GITCOMMIT := $(shell git rev-parse --short HEAD)
TAGS := $(shell git tag --contains $(GITCOMMIT))
DTAG = $(shell git tag -l [0-9][0-9][0-1][0-9][0-3][0-9][0-2][0-9][0-5][0-9]-* | tail -n1 | sed 's/^[0-9]\{10\}-//')


TOP_DIR = $(shell python -c "import os.path as p; print p.abspath('../..')")

TOP_DIR_NAME = $(shell basename $(TOP_DIR))

DIR = $(shell pwd)

LIB_DIR = lib

LBIN_DIR = bin

EXECUTABLE_SCRIPT_NAME = run_$(MODULE_CAPS).sh

default: compile build-executable-script-python

compile:
	kb-sdk compile $(SPEC_FILE) \
		--out $(LIB_DIR) \
		--pyclname biokbase.$(MODULE).$(MODULE_CAPS)Client \
		--pysrvname biokbase.$(MODULE).$(MODULE_CAPS) \
		--pyimplname biokbase.$(MODULE).$(MODULE_CAPS)Impl;

# NOTE: script generation and wrapping in various languages should be
# handled in a kb-mobu tool, but for now we just generate the
# script within this makefile
build-executable-script-python: setup-local-dev-kb-py-libs
	mkdir -p $(LBIN_DIR)
	echo '#!/bin/bash' > $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_DEPLOYMENT_CONFIG="$(DIR)/deploy.cfg"' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_SERVICE_NAME="$(MODULE_CAPS)"' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export PYTHONPATH="$(DIR)/$(LIB_DIR)"' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'python $(DIR)/lib/biokbase/$(MODULE)/$(MODULE_CAPS).py $$1 $$2 $$3' \
		>> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	chmod +x $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
ifeq ($(TOP_DIR_NAME), dev_container)
	cp $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME) $(TOP_DIR)/bin/.
else
	cp $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME) $(KB_TOP)/bin/.
endif


##
# WARNING Hardcoded /kb/dev_container (shouldn't be this way...)
#
setup-local-dev-kb-py-libs:
	touch lib/biokbase/__init__.py
	touch lib/biokbase/$(MODULE)/__init__.py
	rsync -vrh /kb/dev_container/modules/kbapi_common/lib/biokbase/* lib/biokbase/.
	rsync -vrh /kb/dev_container/modules/auth/lib/biokbase/* lib/biokbase/.
	rsync -vrh /kb/dev_container/modules/handle_service/lib/biokbase/* lib/biokbase/.
	rsync -vrh /kb/dev_container/modules/workspace_deluxe/lib/biokbase/* lib/biokbase/.

clean:
	rm -rfv $(LBIN_DIR)



# below are targets for deploying in a KBase environment - note that these
# are hacked together to get things working for now, and should be refactored if
# this example is going to be copied into a production service
ifeq ($(TOP_DIR_NAME), dev_container)
include $(TOP_DIR)/tools/Makefile.common
include $(TOP_DIR)/tools/Makefile.common.rules

DEPLOY_RUNTIME ?= /kb/runtime
TARGET ?= /kb/deployment
#SERVICE_DIR ?= $(TARGET)/services/$(MODULE)

deploy: deploy-scripts

deploy-scripts: deploy-libs2 deploy-executable-script
	bash $(DIR)/deps/pylib.sh

deploy-service: deploy-libs2 deploy-executable-script deploy-service-scripts deploy-cfg

deploy-libs2:
	kb-sdk install AssemblyUtil
	@echo "Deploying libs to target: $(TARGET)"
	mkdir -p $(TARGET)/lib/biokbase
	rsync -vrh lib/biokbase/$(MODULE) $(TARGET)/lib/biokbase/.
	rsync -vrh lib/AssemblyUtil $(TARGET)/lib/.
	rsync -vrh lib/GenomeFileUtil $(TARGET)/lib/.

deploy-executable-script:
	@echo "Installing executable scripts to target: $(TARGET)/bin"
	echo '#!/bin/bash' > $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_RUNTIME=$(DEPLOY_RUNTIME)' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export PYTHONPATH="$(TARGET)/lib"' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'python $(TARGET)/lib/biokbase/$(MODULE)/$(MODULE_CAPS).py $$1 $$2 $$3' \
		>> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	chmod +x $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)

deploy-service-scripts:
	@echo "Preparing start/stop scripts for service"

deploy-cfg:
	@echo "Generating real deployment.cfg based on template"

test: test-impl create-test-wrapper


test-impl: create-test-wrapper
	./test/script_test/run_tests.sh

create-test-wrapper:
	@echo "Creating test script wrapper in test/script_test"
	echo '#!/bin/bash' > test/script_test/run_tests.sh
	echo 'export KB_RUNTIME=$(DEPLOY_RUNTIME)' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME) >> test/script_test/run_tests.sh
	echo 'export PYTHONPATH="$(DIR)/$(LIB_DIR)"' >> test/script_test/run_tests.sh
	echo 'export KB_DEPLOYMENT_CONFIG="$(DIR)/deploy.cfg"' >> test/script_test/run_tests.sh
	echo 'python $(DIR)/test/script_test/basic_test.py $$1 $$2 $$3' \
		>> test/script_test/run_tests.sh
	chmod +x test/script_test/run_tests.sh

else
####
# Assumption: make deploy is done before test, which will be satisfied in Dockerfile
####
# TEMP HACK
KB_TOP = /kb/dev_container
TOP_DIR = $(KB_TOP)
KB_SERVICE_NAME = $(MODULE_CAPS)
KB_DEPLOYMENT_CONFIG = $(TARGET)/deployment.cfg
include $(KB_TOP)/tools/Makefile.common
include $(KB_TOP)/tools/Makefile.common.rules

DEPLOY_RUNTIME ?= /kb/runtime
TARGET ?= /kb/deployment
#SERVICE_DIR ?= $(TARGET)/services/$(MODULE)

deploy: deploy-lscripts deploy-cfg
	

deploy-ensure-dirs:
	if [ ! -d $(TARGET)/shbin ]; then rm -rf $(TARGET)/shbin; mkdir -p $(TARGET)/shbin; fi

deploy-lscripts: deploy-ensure-dirs deploy-libs2 deploy-executable-script deploy-scripts
	bash $(DIR)/deps/pylib.sh

deploy-service: deploy-libs2 deploy-executable-script deploy-service-scripts deploy-cfg

deploy-libs2:
	@echo "Deploying libs to target: $(TARGET)"
	mkdir -p $(TARGET)/lib/biokbase
	rsync -vrh lib/biokbase/$(MODULE) $(TARGET)/lib/biokbase/.
	rsync -vrh lib/AssemblyUtil $(TARGET)/lib/.
	rsync -vrh lib/GenomeFileUtil $(TARGET)/lib/.

deploy-executable-script:
	@echo "Installing executable scripts to target: $(TARGET)/bin"
	echo '#!/bin/bash' > $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_RUNTIME=$(DEPLOY_RUNTIME)' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export PYTHONPATH="$(TARGET)/lib"' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_SERVICE_NAME="$(MODULE_CAPS)"' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export KB_DEPLOYMENT_CONFIG="$(DIR)/deploy.cfg"' >> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	echo 'python $(TARGET)/lib/biokbase/$(MODULE)/$(MODULE_CAPS).py $$1 $$2 $$3' \
		>> $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)
	chmod +x $(TARGET)/bin/$(EXECUTABLE_SCRIPT_NAME)

deploy-service-scripts:
	@echo "TODO: Preparing start/stop scripts for service"

#deploy-cfg:
#	@echo "TODO: Generating real deployment.cfg based on template"

#test: deploy test-impl create-test-wrapper 
test: test-impl create-test-wrapper

test-impl: create-test-wrapper
	./test/script_test/run_tests.sh

create-test-wrapper:
	@echo "Creating test script wrapper in test/script_test"
	echo '#!/bin/bash' > test/script_test/run_tests.sh
	echo 'export KB_RUNTIME=$(DEPLOY_RUNTIME)' >> test/script_test/run_tests.sh
	echo 'export PYTHONPATH="$(TARGET)/lib"' >> test/script_test/run_tests.sh
	echo 'export KB_SERVICE_NAME="$(MODULE_CAPS)"' >> test/script_test/run_tests.sh
	echo 'export KB_DEPLOYMENT_CONFIG="$(DIR)/deploy.cfg"' >> test/script_test/run_tests.sh
	echo 'python $(DIR)/test/script_test/basic_test.py $$1 $$2 $$3' \
		>> test/script_test/run_tests.sh
	chmod +x test/script_test/run_tests.sh

endif

