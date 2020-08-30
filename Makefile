all:
	$(MAKE) -C src

lib:
	$(MAKE) -C src/galib/src

clean:
	$(MAKE) -C src clean
