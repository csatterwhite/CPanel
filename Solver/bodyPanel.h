//
//  bodyPanel.h
//  CPanel
//
//  Created by Chris Satterwhite on 5/1/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__bodyPanel__
#define __CPanel__bodyPanel__

#include <iostream>
#include "panel.h"

class bodyPanel : public panel
{
    double sigma;
    double mu;
    
public:
    bodyPanel(short surfaceID) : panel(surfaceID) {}
};

#endif /* defined(__CPanel__bodyPanel__) */
