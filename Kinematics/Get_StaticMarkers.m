function [staticR, staticMarker] = Get_StaticMarkers(Markers)

[staticR, ASISBreadth, HipOrigin, staticJC] = HipRotationMatrix(Markers);

[staticRMK, localRKNE, RThigh_R, staticJC] = nonAnan_RMatrices('RThigh', Markers, staticJC);
staticR.RThigh = RThigh_R;
staticMarker.RMK = staticRMK;
staticMarker.RKNE = localRKNE;

[staticRMA, localRANK, RShank_R, staticJC] = nonAnan_RMatrices('RShank', Markers, staticJC);
staticR.RShank = RShank_R;
staticMarker.RMA = staticRMA;
staticMarker.RANK = localRANK;

[staticLMK, localLKNE, LThigh_R, staticJC] = nonAnan_RMatrices('LThigh', Markers, staticJC);
staticR.LThigh = LThigh_R;
staticMarker.LMK = staticLMK;
staticMarker.LKNE = localLKNE;

[staticLMA, localLANK, LShank_R, staticJC] = nonAnan_RMatrices('LShank', Markers, staticJC);
staticR.LShank = LShank_R;
staticMarker.LMA = staticLMA;
staticMarker.LANK = localLANK;

end

