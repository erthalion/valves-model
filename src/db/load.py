#!/usr/bin/python

from pymongo import MongoClient
import gridfs
import ConfigParser
import sys

config = ConfigParser.RawConfigParser()
config.read('mongolab.config')

def load_files(result, area, description, database):
    """ Load the result data file and associated area data file into the database
    Args:
        result (type): result of the calculations
        area (type): it's described calculation area
        database (type): destination
    """
    gridFs = gridfs.GridFS(database)
    area_id = gridFs.put( area.read(), filename="area", description=description)
    result_id = gridFs.put( result.read(), filename="result", area=area_id, description=description)

def main():
    if len(sys.argv) < 2:
        sys.exit('Usage: %s description' % sys.argv[0])

    area = open("mask.vtk", 'r')
    result = open("streamlines.vtk", 'r')

    connection = MongoClient(config.get('database', 'connection_string'))
    database = connection.valves
    load_files(result, area, sys.argv[1], database)

if __name__ == '__main__':
    main()
