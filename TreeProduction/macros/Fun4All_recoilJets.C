//======================================================================
//  Fun4All_recoilJets.C
//  --------------------------------------------------------------------
#pragma once
#define RJ_UNIFIED_ANALYSIS_PP 1
#include "Fun4All_recoilJets_unified_impl.C"

void Fun4All_recoilJets(const int   nEvents   =  0,
                        const char* listFile  = "input_files.list",
                        const char* outRoot   = "TrigPlot.root",
                        const bool  verbose   = false)
{
  Fun4All_recoilJets_unified_impl(nEvents, listFile, outRoot, verbose);
}
