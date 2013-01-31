#!/usr/bin/python

import pymongo
import gridfs

def load_files(result, area, database):
    """ Load the result data file and associated area data file into the database
    Args:
        result (type): result of the calculations
        area (type): it's described calculation area
        database (type): destination
    """
    gridFs = gridfs.GridFS(database)
    area_id = gridFs.put( area.read(), filename="area")
    result_id = gridFs.put( result.read(), filename="result", area=area_id)

def main():
    area = open("mask.vtk", 'r')
    result = open("streamlines.vtk", 'r')

    connection = pymongo.Connection( "localhost", 27017)
    database = connection.valves
    load_files(result, area, database)

if __name__ == '__main__':
    main()
