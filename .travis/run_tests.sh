#!/bin/bash

coverage run --source=sbms --omit="sbms/tests" setup.py test
