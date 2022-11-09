IMG_DIR=mkdocs/img
OUT_DIR=site
PYTHON=python
EXAMPLES = $(wildcard ../examples/ex*.py)
EXAMPLE_IMGS = ${EXAMPLES:.py=.png}

all: mkdocs

mkdocs: $(EXAMPLE_IMGS)
	mkdocs build -d $(OUT_DIR)

$(EXAMPLE_IMGS): $(EXAMPLES)
	@mkdir -p $(IMG_DIR)
	@echo "Building example figures..."
	@for py in $^; do $(PYTHON) $$py; done
	@echo "Copying images to output directory..."
	@for img in $(EXAMPLE_IMGS); do cp $$img $(IMG_DIR); done

clean:
	rm -rf $(IMG_DIR)
	rm -rf $(OUT_DIR)
	@for img in $(EXAMPLE_IMGS); do rm -f $$img; done
