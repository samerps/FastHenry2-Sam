# Master makefile for fastcap and related programs and documentation

SRC = src
FASTH = $(SRC)/fasthenry

default:
	@echo Specify what to make:
	@echo " fasthenry - inductance calculation program"
	@echo " clean - remove object files"

fasthenry:
	cd $(FASTH) ; $(MAKE) fasthenry


clean:
	cd $(FASTH) ; $(MAKE) clean


