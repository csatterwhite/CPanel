//
//  member.h
//  CPanel
//
//  Created by Chris Satterwhite on 4/24/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef CPanel_member_h
#define CPanel_member_h

template<typename type>
class member
{
    typedef std::array<double,3>  point;
    type* pObj;
    point refPoint;
    
public:
    member(type* object, point referencePoint) : pObj(object), refPoint(referencePoint) {}
    
    point getRefPoint() const {return refPoint;}
    
    type *getObj() {return pObj;}
};


#endif
