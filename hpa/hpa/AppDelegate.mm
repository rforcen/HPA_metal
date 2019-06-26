//
//  AppDelegate.m
//  hpa
//
//  Created by asd on 24/06/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#import "AppDelegate.h"
#import "hpa/hpa.hpp"
#import "metal/MetalDevice.h"
#import "Thread.h"



@interface AppDelegate ()

@end

@implementation AppDelegate {
    HPA x, y;
}

-(NSTimeInterval) timeIt: (void (^) (void))block {
    NSDate *start = [NSDate date];
    block();
    return -[start timeIntervalSinceNow];
}

-(void) testMetal {
    int __n=400, w=__n, h=__n, n=w*h, size=n*sizeof(HPA);
    HPA*hpasMetal, *hpasST=new HPA[n], *hpasMT=new HPA[n];
    
    printf("metal/cpu HPA comp. for %d items\n", n);
    
    MetalDevice*md=[MetalDevice init];
    hpasMetal=(HPA*)[md newData:size];
    memset(hpasMetal, 0, size);

    [md compileFunc:@"testHPA"];
    id<MTLBuffer> hpasBuff=[md copyBuffer:hpasMetal length:size];
    [md setBufferParam:hpasBuff index:0];
    [md setIntParam:&w index:1];
    
    NSTimeInterval tgpu = [MetalDevice timeIt:^{
        [md runThreadsWidth:w height:h];
    }];
    
    NSTimeInterval tst=[MetalDevice timeIt:^{ // single thread
        for (int i=0; i<n; i++) {
            HPA x=i;
            hpasST[i]=x*x;
        }
    }];
    
    NSTimeInterval tmt=[MetalDevice timeIt:^{ // multi thread
        Thread(n).run([&hpasMT](int t, int from, int to){
            for (int i=from; i<to; i++) {
                HPA x=i;
                hpasMT[i]=x*x;
            }
        });
    }];
    
    // check ST/MT
    HPA s=0;
    for (int i=0; i<n; i++) s+=abs(hpasST[i]-hpasMT[i]);
    printf("difference st/mt: %s -> %s\n", s.tochar(), s==0 ? "OK":" non EQ!");
    
    // check ST/gpu
    s=0;
    for (int i=0; i<n; i++) s+=abs(hpasST[i]-hpasMetal[i]);
    printf("difference st/metal: %s -> %s\n", s.tochar(), s==0 ? "OK":" non EQ!");
    
    // timings 0..10 ,...., n-10..n
    printf("time gpu: %.1f, time st:%.1f, mt:%.1f, ratio st/mt: %.2f\n\n", tgpu, tst, tmt, tst/tmt);
    for (int i=0; i<10; i++) printf("%s\t", hpasMetal[i].tochar(6));
    puts("\n......");
    for (int i=n-10; i<n; i++) printf("%s\t", hpasMetal[i].tochar(6));
    puts("\n\nend");
    
    // release mem
    [md deleteData:size data:hpasMetal];
    delete[]hpasST;
    delete[]hpasMT;
}

-(void) testHPA {
    HPA x0(5.5), x1(2.3);
    printf("%s+%s=%s, *=%s, /=%s\n", x1.tochar(), x0.tochar(), (x1+x0).tochar(), (x0*x1).tochar(), (-x0/x1).tochar(30));
    
    
    x = 456;
    y = "1.23678111222333444555666777888999e4230";
    
    printf("y=%s\n", y.tochar(39));
    
    NSTimeInterval thpa = [self timeIt:^{
        
        HPA y("123.345e1234"), z, s;
        
        for (int i=0; i<10000; i++) {
            z+=y;
            self->x-=z;
            s+=sin(self->x) * cos(self->x) + tan(self->x);
        }
        
        char*sy=self->x.tochar();
        printf("x=%s y=%s z=%s, s=%s, sy=%s",self->x.tochar(), y.tochar(), z.tochar(), s.tochar(), sy);
    }];
    
    printf("\nlap time: %f\n", thpa);
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    [self testMetal];
}

- (BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)sender { return YES;}
- (void)applicationWillTerminate:(NSNotification *)aNotification {}


@end
