//======================================================================
//  Fun4All_recoilJets_AuAu.C
//  --------------------------------------------------------------------
#pragma once
#define RJ_UNIFIED_ANALYSIS_AUAU 1
#include "Fun4All_recoilJets_unified_impl.C"

void Fun4All_recoilJets_AuAu(const int   nEvents   =  0,
                             const char* listFile  = "input_files.list",
                             const char* outRoot   = "TrigPlot.root",
                             const bool  verbose   = false)
{
  Fun4All_recoilJets_unified_impl(nEvents, listFile, outRoot, verbose);
}
