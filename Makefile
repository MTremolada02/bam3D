CXX = g++
GFALIBS_DIR := $(CURDIR)/gfalibs
INCLUDE_DIR = -I./include -I$(GFALIBS_DIR)/include #-Imodule2/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS) $(CFLAGS)

TARGET = bam3D #name of tool
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o
LIBS = -lz -lhts
LDFLAGS = -pthread

OBJS := main runner
BINS := $(addprefix $(BINDIR)/, $(OBJS))

head: $(BINS) bam3D | $(BUILD)#module2 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(GFALIBS_DIR)/*.o $(LIBS)
	
debug: CXXFLAGS += -DDEBUG -O0
debug: CCFLAGS += -DDEBUG
debug: head

all: head validate regenerate

$(OBJS): %: $(BINDIR)/%
	@
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@
	
.PHONY: bam3D
bam3D:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"

validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/bam3D-validate $(SOURCE)/$(TEST_TARGET).cpp

regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/bam3D-generate-tests $(SOURCE)/$(GENERATE_TARGET).cpp

#.PHONY: module2
#module2:
#	$(MAKE) -j -C $(MODULE2_DIR) CXXFLAGS="$(CXXFLAGS)"
	
$(BUILD):
$(BINDIR):
	-mkdir -p $@

clean:
	$(MAKE) -j -C $(GFALIBS_DIR) clean
#	$(MAKE) -j -C $(MODULE2_DIR) clean
	$(RM) -r build
