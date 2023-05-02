all:
	@gcc -o sssp_serial src/util.h src/util.c  src/sssp_serial.c 
#	./sssp_serial input/sample_graph_paper_csc_8v_16e_weighted.txt 6 out.txt

clear:
	@rm -f a.out out.txt sssp_serial

c: clear

z: zip

zip:
	@zip -r muhammed_can_kucukaslan.zip Makefile src/
