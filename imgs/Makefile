SVGFILES := $(wildcard *.svg)

all: $(SVGFILES:%.svg=%.pdf)

%.pdf : %.svg
	inkscape -D -z --export-pdf=$(CURDIR)/$@ $(CURDIR)/$<
