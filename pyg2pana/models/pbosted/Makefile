MODULE := pbosted
EXT_SUFFIX := $(shell python3-config --extension-suffix)

all: _$(MODULE)$(EXT_SUFFIX)

_$(MODULE)$(EXT_SUFFIX): pbosted.f95 F1F209.f pbosted.pyf
	f2py -c pbosted.pyf pbosted.f95 F1F209.f

clean:
	rm -rf _$(MODULE)$(EXT_SUFFIX)
	rm -rf _$(MODULE)$(EXT_SUFFIX).*