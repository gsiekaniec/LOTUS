#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

def true_stem(path):
	'''
	Take a file Path and return the true stem of this path (file name without extensions)
	'''
	stem = Path(path).stem
	return stem if stem == path else true_stem(stem)
