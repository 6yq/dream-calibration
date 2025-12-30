APP:=dream
ENERGIES:=20 40 80 100 150 200 250 300
N_EVENTS:=1000
PARTICLE:=pi-

SRC_DIR:=src
INC_DIR:=include
ROOT_MAC:=root/analysis.c

BUILD_DIR:=build
OBJ_DIR:=$(BUILD_DIR)/obj
BIN:=$(BUILD_DIR)/$(APP)

SIM_DIR:=sim
MAC_DIR:=mac
RES_DIR:=results

VIS_MAC:=visVRML.mac

CXX:=g++
CXXFLAGS:=-O2 -std=c++17 -I$(INC_DIR) $(shell geant4-config --cflags) $(shell root-config --cflags)
LDLIBS:=$(shell geant4-config --libs) $(shell root-config --libs)

SRCS:=$(wildcard $(SRC_DIR)/*.cc)
OBJS:=$(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SRCS))

.PHONY: all builds vis sims analyze clean

all: analyze

$(BUILD_DIR) $(OBJ_DIR) $(SIM_DIR) $(MAC_DIR) $(RES_DIR):
	mkdir -p $@

builds: $(BIN)
$(BIN): $(OBJS) | $(BUILD_DIR)
	$(CXX) $^ -o $@ $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

vis: builds
	$(BIN) -m $(VIS_MAC) > $@.log

SIM_OUT:=$(addprefix $(SIM_DIR)/E,$(addsuffix .root,$(ENERGIES)))
sims: $(SIM_OUT)

$(MAC_DIR)/run_E%.mac: | $(MAC_DIR) $(SIM_DIR)
	echo "/run/initialize" > $@
	echo "/control/verbose 1" >> $@
	echo "/run/verbose 1" >> $@
	echo "/event/verbose 0" >> $@
	echo "/tracking/verbose 0" >> $@
	echo "/gun/particle $(PARTICLE)" >> $@
	echo "/gun/energy $* GeV" >> $@
	echo "/ana/output $(SIM_DIR)/E$*.root" >> $@
	echo "/run/beamOn $(N_EVENTS)" >> $@

$(SIM_DIR)/E%.root: $(BIN) $(MAC_DIR)/run_E%.mac | $(SIM_DIR)
	$(BIN) -m $(MAC_DIR)/run_E$*.mac > $@.log

analyze: $(SIM_OUT) | $(RES_DIR)
	root -l -b -q '$(ROOT_MAC)("$(SIM_DIR)","$(RES_DIR)")'

clean:
	rm -rf $(BUILD_DIR) $(SIM_DIR) $(RES_DIR)
	rm -f $(MAC_DIR)/run_E*.mac
