#COMPILER MODE C++11
CXX=g++ -std=c++11

#COMPILER FLAGS
CXXFLAG_REL=-O2
CXXFLAG_DBG=-g
CXXFLAG_WRN=-Wall -Wextra -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated -Wno-unused-parameter

#LINKER FLAGS
LDFLAG_REL=-O2
LDFLAG_DBG=-g

#BASE LIBRARIES
LIB_FLAGS=-lm -lz -lgsl -lblas -lbz2 -lpthread

#FILE LISTS
BFILE=bin/clomics
HFILE=$(shell find src -name *.h)
TFILE=$(shell find lib -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#DEFAULT VERSION (I.E. UNIGE DESKTOP RELEASE VERSION)
all: desktop

#UNIGE DESKTOP RELEASE VERSION
desktop: RMATH_INC=/home/olivier/Tools/R-3.2.2/src/include
desktop: RMATH_LIB=/home/olivier/Tools/R-3.2.2/src/nmath/standalone
desktop: HTSLD_INC=/home/olivier/Tools/htslib-1.3
desktop: HTSLD_LIB=/home/olivier/Tools/htslib-1.3
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB=/usr/lib/x86_64-linux-gnu
desktop: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
desktop: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
desktop: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
desktop: LDFLAG=$(LDFLAG_REL)
desktop: $(BFILE)

#UNIGE DESKTOP DEBUG VERSION
desktop-dbg: RMATH_INC=/home/olivier/Tools/R-3.2.2/src/include
desktop-dbg: RMATH_LIB=/home/olivier/Tools/R-3.2.2/src/nmath/standalone
desktop-dbg: HTSLD_INC=/home/olivier/Tools/htslib-1.3
desktop-dbg: HTSLD_LIB=/home/olivier/Tools/htslib-1.3
desktop-dbg: BOOST_INC=/usr/include
desktop-dbg: BOOST_LIB=/usr/lib/x86_64-linux-gnu
desktop-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
desktop-dbg: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
desktop-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
desktop-dbg: LDFLAG=$(LDFLAG_DBG)
desktop-dbg: $(BFILE)

#UNIGE SERVER RELEASE VERSION (=BINARY RELEASED TO PUBLIC)
server: RMATH_INC=/home/popgen/delaneau/SOFT/R-3.2.1/src/include
server: RMATH_LIB=/home/popgen/delaneau/SOFT/R-3.2.1/src/nmath/standalone
server: HTSLD_INC=/home/popgen/delaneau/SOFT/htslib-1.2.1
server: HTSLD_LIB=/home/popgen/delaneau/SOFT/htslib-1.2.1
server: BOOST_INC=/home/popgen/delaneau/SOFT/boost/boost_1_59_0/build/include
server: BOOST_LIB=/home/popgen/delaneau/SOFT/boost/boost_1_59_0/build/lib
server: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
server: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
server: LIB_FLAGS=$(LIB_FLAGS) -lgslcblas
server: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
server: LDFLAG=$(LDFLAG_REL)
server: $(BFILE)

#DELL LAPTOP RELEASE VERSION
laptop: RMATH_INC=/home/olivier/Libraries/R-3.2.2/src/include
laptop: RMATH_LIB=/home/olivier/Libraries/R-3.2.2/src/nmath/standalone
laptop: HTSLD_INC=/home/olivier/Libraries/htslib-1.2.1
laptop: HTSLD_LIB=/home/olivier/Libraries/htslib-1.2.1
laptop: BOOST_INC=/usr/include
laptop: BOOST_LIB=/usr/lib/x86_64-linux-gnu
laptop: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
laptop: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
laptop: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
laptop: LDFLAG=$(LDFLAG_REL)
laptop: $(BFILE)

#DELL LAPTOP DEBUG VERSION
laptop-dbg: RMATH_INC=/home/olivier/Libraries/R-3.2.2/src/include
laptop-dbg: RMATH_LIB=/home/olivier/Libraries/R-3.2.2/src/nmath/standalone
laptop-dbg: HTSLD_INC=/home/olivier/Libraries/htslib-1.2.1
laptop-dbg: HTSLD_LIB=/home/olivier/Libraries/htslib-1.2.1
laptop-dbg: BOOST_INC=/usr/include
laptop-dbg: BOOST_LIB=/usr/lib/x86_64-linux-gnu
laptop-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
laptop-dbg: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
laptop-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
laptop-dbg: LDFLAG=$(LDFLAG_DBG)
laptop-dbg: $(BFILE)

#VITAL-IT RELEASE VERSION
cluster: RMATH_INC=/software/R/3.1.1/include
cluster: RMATH_LIB=/software/R/3.1.1/lib64
cluster: HTSLD_INC=/software/UHTS/Analysis/samtools/1.2/include
cluster: HTSLD_LIB=/software/UHTS/Analysis/samtools/1.2/lib64
cluster: BOOST_INC=/software/include
cluster: BOOST_LIB=/software/lib64
cluster: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
cluster: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
cluster: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
cluster: LDFLAG=$(LDFLAG_REL)
cluster: $(BFILE)

#VITAL-IT DEBUG VERSION
cluster-dbg: RMATH_INC=/software/R/3.1.1/include
cluster-dbg: RMATH_LIB=/software/R/3.1.1/lib64
cluster-dbg: HTSLD_INC=/software/UHTS/Analysis/samtools/1.2/include
cluster-dbg: HTSLD_LIB=/software/UHTS/Analysis/samtools/1.2/lib64
cluster-dbg: BOOST_INC=/software/include
cluster-dbg: BOOST_LIB=/software/lib64
cluster-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
cluster-dbg: IFLAG=-Ilib/OTools -Ilib -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
cluster-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
cluster-dbg: LDFLAG=$(LDFLAG_DBG)
cluster-dbg: $(BFILE)

#UBUNTU RELEASE VERSION
ubuntu: RMATH_INC=$(HOME)/R-3.2.2/src/include
ubuntu: RMATH_LIB=$(HOME)/R-3.2.2/src/nmath/standalone
#ubuntu: RMATH_INC=/usr/include
#ubuntu: RMATH_LIB=/usr/lib
ubuntu: HTSLD_INC=/usr/local/include/
ubuntu: HTSLD_LIB=/usr/local/lib
ubuntu: BOOST_INC=/usr/include
ubuntu: BOOST_LIB=/usr/lib/x86_64-linux-gnu
ubuntu: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
ubuntu: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
ubuntu: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
ubuntu: LDFLAG=$(LDFLAG_REL)
ubuntu: $(BFILE)

#UBUNTU DEBUG VERSION
ubuntu-dbg: RMATH_INC=$(HOME)/R-3.2.2/src/include
ubuntu-dbg: RMATH_LIB=$(HOME)/R-3.2.2/src/nmath/standalone
#ubuntu-dbg: RMATH_INC=/usr/include
#ubuntu-dbg: RMATH_LIB=/usr/lib
ubuntu-dbg: HTSLD_INC=/usr/local/include/
ubuntu-dbg: HTSLD_LIB=/usr/local/lib
ubuntu-dbg: BOOST_INC=/usr/include
ubuntu-dbg: BOOST_LIB=/usr/lib/x86_64-linux-gnu
ubuntu-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
ubuntu-dbg: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
ubuntu-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a
ubuntu-dbg: LDFLAG=$(LDFLAG_DBG)
ubuntu-dbg: $(BFILE)

#MAC RELEASE VERSION
mac: RMATH_INC=$(HOME)/Libraries/R-3.2.2/src/include
mac: RMATH_LIB=$(HOME)/Libraries/R-3.2.2/src/nmath/standalone
mac: HTSLD_INC=$(HOME)/Libraries/htslib-1.2.1
mac: HTSLD_LIB=$(HOME)/Libraries/htslib-1.2.1
mac: BOOST_INC=/opt/local/include
mac: BOOST_LIB=/opt/local/lib
mac: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
mac: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
mac: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams-mt.a $(BOOST_LIB)/libboost_program_options-mt.a
mac: LDFLAG=$(LDFLAG_REL) -L /opt/local/lib
mac: $(BFILE)

#MAC DEBUG VERSION
mac-dbg: RMATH_INC=$(HOME)/Libraries/R-3.2.2/src/include
mac-dbg: RMATH_LIB=$(HOME)/Libraries/R-3.2.2/src/nmath/standalone
mac-dbg: HTSLD_INC=$(HOME)/Libraries/htslib-1.2.1
mac-dbg: HTSLD_LIB=$(HOME)/Libraries/htslib-1.2.1
mac-dbg: BOOST_INC=/opt/local/include
mac-dbg: BOOST_LIB=/opt/local/lib
mac-dbg: CXXFLAG=$(CXXFLAG_DBG) $(CXXFLAG_WRN)
mac-dbg: IFLAG=-Ilib/OTools -Ilib/ -I$(RMATH_INC) -I$(HTSLD_INC) -I$(BOOST_INC)
mac-dbg: LIB_FILES=$(RMATH_LIB)/libRmath.a $(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams-mt.a $(BOOST_LIB)/libboost_program_options-mt.a
mac-dbg: LDFLAG=$(LDFLAG_DBG) -L /opt/local/lib
mac-dbg: $(BFILE)


#COMPILATION RULES
$(BFILE): $(OFILE)
	$(CXX) $^ $(LIB_FILES) -o $@ $(LIB_FLAGS) $(LDFLAG)

obj/clomics.o: src/clomics.cpp $(HFILE) $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/data.o: src/common/data.cpp src/common/data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/build_%.o: build_%.cpp build_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/score_%.o: score_%.cpp score_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/locate_%.o: locate_%.cpp locate_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/call_%.o: call_%.cpp call_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/topo_%.o: topo_%.cpp topo_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

obj/corr_%.o: corr_%.cpp corr_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/quantify_%.o: quantify_%.cpp quantify_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/qtl_%.o: qtl_%.cpp qtl_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/call_%.o: call_%.cpp call_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/hic_%.o: hic_%.cpp hic_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/motif_cis_%.o: motif_cis_%.cpp motif_cis_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/motif_trans_%.o: motif_trans_%.cpp motif_trans_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/trans_%.o: trans_%.cpp trans_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/span_%.o: span_%.cpp span_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/trans_filter_%.o: trans_filter_%.cpp trans_filter_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/stat_%.o: stat_%.cpp stat_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/mbed_%.o: mbed_%.cpp mbed_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/merge_%.o: merge_%.cpp merge_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/reptime_%.o: reptime_%.cpp reptime_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/overlap_%.o: overlap_%.cpp overlap_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/gpeak_%.o: gpeak_%.cpp gpeak_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/vcfcollapse_%.o: vcfcollapse_%.cpp vcfcollapse_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)
	
obj/binding_%.o: binding_%.cpp binding_data.h $(TFILE)
	$(CXX) -o $@ -c $< $(CXXFLAG) $(IFLAG)

clean: 
	rm -f obj/*.o $(BFILE)