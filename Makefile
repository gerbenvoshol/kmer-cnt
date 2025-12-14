CFLAGS=-g -Wall -O2
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz
MATHLIBS=-lm
HTSLIB=htslib/libhts.a
HTSINC=-Ihtslib
HTSLIBS=$(HTSLIB) -lz -lm -lpthread
PROG=kc-c1 kc-c2 kc-c3 kc-c4 kc-cpp1 kc-cpp2 yak-count snp-pattern-gen vaf-counter correlation-matrix match-classifier bam-vaf-counter vcf-vaf-counter

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)

kc-c1:kc-c1.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

kc-c2:kc-c2.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

kc-c3:kc-c3.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ kc-c3.c kthread.c $(LIBS) -lpthread

kc-c4:kc-c4.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ kc-c4.c kthread.c $(LIBS) -lpthread

yak-count:yak-count.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ yak-count.c kthread.c $(LIBS) -lpthread

kc-cpp1:kc-cpp1.cpp ketopt.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

kc-cpp2:kc-cpp2.cpp ketopt.h robin_hood.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

snp-pattern-gen:snp-pattern-gen.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

vaf-counter:vaf-counter.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ vaf-counter.c kthread.c $(LIBS) -lpthread

correlation-matrix:correlation-matrix.c ketopt.h
	$(CC) $(CFLAGS) -o $@ $< $(MATHLIBS)

match-classifier:match-classifier.c ketopt.h
	$(CC) $(CFLAGS) -o $@ $< $(MATHLIBS)

$(HTSLIB):
	cd htslib && ./configure --disable-bz2 --disable-lzma && $(MAKE) lib-static

bam-vaf-counter:bam-vaf-counter.c khashl.h ketopt.h kthread.h $(HTSLIB)
	$(CC) $(CFLAGS) $(HTSINC) -o $@ bam-vaf-counter.c kthread.c $(HTSLIBS) -lcurl -lcrypto -ldeflate

vcf-vaf-counter:vcf-vaf-counter.c ketopt.h $(HTSLIB)
	$(CC) $(CFLAGS) $(HTSINC) -o $@ vcf-vaf-counter.c $(HTSLIBS) -lcurl -lcrypto -ldeflate

clean:
	rm -fr *.dSYM $(PROG)
	cd htslib && $(MAKE) clean 2>/dev/null || true
