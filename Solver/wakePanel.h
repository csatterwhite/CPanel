//
//  wakePanel.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wakePanel__
#define __CPanel__wakePanel__

#include <iostream>
#include "panel.h"

class wakePanel : public panel
{
    double mu;
    
public:
    wakePanel(short surfaceID) : panel(surfaceID) {}
};

#endif /* defined(__CPanel__wakePanel__) */
