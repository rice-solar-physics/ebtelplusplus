API_DIR=mkdocs/api
DOXYGEN_DIR=_tmp_doxygen
IMG_DIR=img
PYTHON=python

mkdocs: examples 

examples: img_dir


img_dir:
	mkdir $(IMG_DIR)

api: xml markdown

markdown: api_dir xml
	moxygen --anchors --output $(API_DIR)/api.md $(DOXYGEN_DIR)/xml

xml: api_dir
	doxygen Doxyfile

api_dir:
	mkdir $(API_DIR)

clean:
	rm -r $(API_DIR)
	rm -r $(DOXYGEN_DIR)