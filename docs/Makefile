PYTHON = python
CLDOC = cldoc

SOURCES = $(wildcard ../source/*.cpp)
HEADERS = $(wildcard ../source/*.h)
EXAMPLES = $(wildcard ../examples/ex*.py)
EXAMPLE_IMGS = ${EXAMPLES:.py=.png}

OUT_DIR = html
MERGE_DIR = ext

DOC_FLAGS = --report --merge $(MERGE_DIR) --output $(OUT_DIR)
CXX_FLAGS = -std=c++11 -I/opt/local/include -I/usr/include/malloc

all: generate

generate-static: $(EXAMPLE_IMGS)
	$(CLDOC) generate $(CXX_FLAGS) -- --static $(DOC_FLAGS) $(SOURCES) $(HEADERS)

generate: $(EXAMPLE_IMGS)
	$(CLDOC) generate $(CXX_FLAGS) -- $(DOC_FLAGS) $(SOURCES) $(HEADERS)

$(EXAMPLE_IMGS): $(EXAMPLES)
	@mkdir -p $(OUT_DIR)
	@echo "Building example figures..."
	@for py in $^; do $(PYTHON) $$py; done
	@echo "Copying images to output directory..."
	@for img in $(EXAMPLE_IMGS); do cp $$img $(OUT_DIR); done

clean:
	@echo "Cleaning up docs..."
	@rm -r $(OUT_DIR)/*
	@echo "Cleaning up example images..."
	@for img in $(EXAMPLE_IMGS); do rm $$img; done

serve: generate
	$(CLDOC) serve $(OUT_DIR)
