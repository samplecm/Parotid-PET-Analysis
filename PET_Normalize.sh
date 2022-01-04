#!/usr/bin/env bash 
dicomautomaton_dispatcher -v \
/media/caleb/WDBlue/PET_PSMA/PSMA_Analysis/SG_PETRT/20/PET_160315/*dcm \
-o CopyImages:ImageSelection=last \
-o ContourWholeImages:ROILabel=everything:ImageSelection=last \
-o NormalizePixels:ImageSelection=last:NormalizedROILabelRegex=.*:ROILabelRegex=.*:Inclusivity=center:ContourOverlap=ignore:Channel=0:Method=suvbw \
-o ThresholdImages:Lower=1:Low=1:Upper=inf:High=inf:Channel=0:ImageSelection=last \
-o SDL_Viewer 