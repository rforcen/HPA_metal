//
//  hpa.metal
//  hpa
//
//  Created by asd on 25/06/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#include <metal_stdlib>
using namespace metal;

#include "hpa_metal.h"

void testAssign(device HPA*hpas) {
    HPA h0=hpas[0], h1=hpas[1];
    h0 += h1;
    assign(hpas[0], h0); // hpas[0]=h0;
}

kernel void testHPA(    device HPA*hpas[[buffer(0)]],
                        const device int&w [[buffer(1)]],
                        uint2 position [[thread_position_in_grid]] )
{
    uint index=position.x + position.y*w;
    HPA x=index;
    assign(hpas[index], x*x); //  return index^2
}
