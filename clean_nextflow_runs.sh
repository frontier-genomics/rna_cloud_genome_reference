#!/bin/bash

# dont delete the temp directory, it holds the batches which take a long time to download
rm -rf .nextflow
rm -rf work/*
rm -f output/*
rm -f temp/*
rm -rf .nextflow.log
rm -rf .nextflow.log.*
