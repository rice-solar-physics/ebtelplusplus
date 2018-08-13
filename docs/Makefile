API_DIR=mkdocs/api
XML_DIR=.
IMG_DIR=mkdocs/img
OUT_DIR=site
PYTHON=python
CLDOC=cldoc
CXX_FLAGS= -std=c++11 -I/opt/local/include -I/usr/include/malloc
DOC_FLAGS=--type xml --output $(XML_DIR)
SOURCES = $(wildcard ../source/*.cpp)
HEADERS = $(wildcard ../source/*.h)
EXAMPLES = $(wildcard ../examples/ex*.py)
EXAMPLE_IMGS = ${EXAMPLES:.py=.png}

all: mkdocs

mkdocs: $(EXAMPLE_IMGS) api
	mkdocs build -d $(OUT_DIR)

$(EXAMPLE_IMGS): $(EXAMPLES)
	@mkdir -p $(IMG_DIR)
	@echo "Building example figures..."
	@for py in $^; do $(PYTHON) $$py; done
	@echo "Copying images to output directory..."
	@for img in $(EXAMPLE_IMGS); do cp $$img $(IMG_DIR); done

api: xml markdown

api_dir:
	mkdir -p $(API_DIR)

xml: api_dir
	$(CLDOC) generate $(CXX_FLAGS) -- $(DOC_FLAGS) $(SOURCES) $(HEADERS)

markdown: api_dir xml
	$(PYTHON) xml_to_md.py --xml_dir $(XML_DIR)/xml --out_dir $(API_DIR)

clean:
	rm -r $(API_DIR)
	rm -r $(XML_DIR)/xml
	rm -r $(IMG_DIR)
	rm -r $(OUT_DIR)
	@for img in $(EXAMPLE_IMGS); do rm $$img; done
