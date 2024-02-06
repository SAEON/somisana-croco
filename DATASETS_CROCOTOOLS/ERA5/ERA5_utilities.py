#!/usr/bin/env python

# Defined utilities to script ERA_request.py 
#
#  Copyright (c) DDONOSO February 2021
#  e-mail:ddonoso@dgeo.udec.cl  


# -------------------------------------------------
# Getting libraries
# -------------------------------------------------
import datetime
import calendar


# -------------------------------------------------
# Increment date by custom months
# -------------------------------------------------
def addmonths4date(date,addmonths):
    month = date.month - 1 + addmonths
    year = date.year + month // 12
    month = month % 12 + 1
    day = min(date.day,calendar.monthrange(year,month)[1])
    return datetime.date(year,month,day)



