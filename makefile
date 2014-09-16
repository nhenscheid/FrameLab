.PHONY: README.md

README.md:
	cd ./src/util/;\
	python ./updateReadme.py;\
	cd ../../;