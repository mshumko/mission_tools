# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:12:10 2017

@author: mike
"""
import numpy as np

def locateConsecutiveNumbers(x):
    """
    NAME:    locateConsecutiveNumbers(x)
    USE:     Locates consecutive numbers in an array.
    INPUT:   An sorted array of numbers.
    RETURNS: startIndex and endIndex arrays of consecutive numbers.
    AUTHOR:  Mykhaylo Shumko
    MOD:     2017-02-15
    """
    conv = np.convolve([1, -1], x, mode = 'valid') - 1
    consecutiveFlag = np.where(conv != 0)[0] + 1
    startInd = np.insert(consecutiveFlag, 0, 0)
    endInd = np.insert(consecutiveFlag, len(consecutiveFlag), len(x)-1)
    return startInd, endInd