//
//  hpa.cpp
//  hpa
//
//  Created by asd on 24/06/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#include "hpa.hpp"

HPA sin(HPA x) { return HPA(x).sin(); }
HPA cos(HPA x) { return HPA(x).cos(); }
HPA tan(HPA x) { return HPA(x).tan(); }
HPA exp(HPA x) { return HPA(x).exp(); }
HPA exp10(HPA x) { return HPA(x).exp10(); }
HPA exp2(HPA x)  { return HPA(x).exp2(); }
HPA log(HPA x)   { return HPA(x).log(); }
HPA log10(HPA x) { return HPA(x).log10(); }
HPA log2(HPA x)  { return HPA(x).log2(); }
HPA pow(HPA x, HPA y) { return HPA(x).pow(y); }
HPA pwr(HPA x, int n) { return HPA(x).pwr(n); }
HPA fmod(HPA x, HPA y) { return HPA(x).fmod(y); }
bool isNaN(HPA x) { return HPA(x).isNaN(); }
bool isInf(HPA x) { return HPA(x).isInf();}
HPA abs(HPA x) { return HPA(x).abs();}
