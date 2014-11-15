//
//  wakeLine.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wakeLine__
#define __CPanel__wakeLine__

#include <stdio.h>
#include "bodyPanel.h"

class wakeLine
{
    bodyPanel* uPan;
    bodyPanel* lPan;
    double y;
    
public:
    wakeLine(bodyPanel* upper, bodyPanel* lower) : uPan(upper), lPan(lower)
    {
        y = (uPan->getCenter()(1)+lPan->getCenter()(1))/2;
    }
    
    double getY() {return y;}
    bodyPanel* getUpper() {return uPan;}
    bodyPanel* getLower() {return lPan;}
};

#endif /* defined(__CPanel__wakeLine__) */
