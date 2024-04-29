import re
import numpy as np
import pandas as pd  # type: ignore

TS_FORMAT = "%Y %m %d %H %M %S"  # 1915 1 13 0 0 0


# THis is taken from
# https://github.com/paknorton/pyPRMS/blob/development/pyPRMS/Streamflow.py
# until pyPRMS is compatiable with pandas > 2.0


class PRMSStreamflowData(object):
    """Class for working with observed streamflow in the PRMS ASCII data file format"""

    def __init__(
        self, filename, missing=-999.0, verbose=False, include_metadata=True
    ):
        self.__missing = missing
        self.filename = filename
        self.__verbose = verbose
        self.__include_metadata = include_metadata

        self.__timecols = 6  # number columns for time in the file
        self.__headercount = None
        self.__metaheader = None
        self.__types = None
        self.__units = {}
        self.__stations = []
        self.__stationIndex = {}  # Lookup of station id to header info
        self.__rawdata = None
        self.__selectedStations = None
        self.__isloaded = False

        self.load_file(self.filename)

    @property
    def data(self):
        """Pandas dataframe of the observed streamflow for each POI"""

        if self.__selectedStations is None:
            return self.__rawdata
        else:
            return self.__rawdata.ix[:, self.__selectedStations]

    @property
    def headercount(self):
        """Number of rows to skip before data begins"""

        return self.__headercount

    @property
    def metaheader(self):
        """List of columns in the metadata section of the data file"""

        return self.__metaheader

    @property
    def numdays(self):
        """The number of days in the period of record"""

        return self.data.shape[0]

    @property
    def size(self):
        """Number of streamgages"""
        return len(self.stations)

    @property
    def stations(self):
        """List of streamgage IDs from streamflow data file"""

        return self.__stations

    @property
    def units(self):
        """Dictionary of units for observation values"""
        return self.__units

    def load_file(self, filename):
        """Read the PRMS ASCII streamflow data file"""

        infile = open(filename, "r")
        rawdata = infile.read().splitlines()
        infile.close()

        it = iter(rawdata)

        self.__headercount = 0

        # We assume if 'ID' and 'Type' header names exist then we have a valid
        # meta-header.
        for line in it:
            self.__headercount += 1

            # Skip lines until we hit the following
            if line[0:10] == "// Station":
                # Read the next line - these are the fieldnames for the station information
                self.__headercount += 1
                self.__metaheader = re.findall(r"[\w]+", next(it))
                break

        cnt = 0
        order = 0  # defines the order of the data types in the dataset
        st = 0

        # Read the station IDs and optional additional metadata
        for line in it:
            self.__headercount += 1

            if line[0:10] == "//////////":
                break

            # Read station information
            # Include question mark in regex below as a valid character since the obs
            # file uses it for missing data in the station information.
            words = re.findall(r"[\w.-]+|[?]", line)  # Break the row up
            curr_fcnt = len(words)

            # Check that number of station information fields remains constant
            if curr_fcnt != len(self.__metaheader):
                if self.__verbose:
                    print(
                        "WARNING: number of header fields changed from %d to %d"
                        % (len(self.__metaheader), curr_fcnt)
                    ),
                    print("\t", words)
                    # exit()

            try:
                if words[self.__metaheader.index("Type")] not in self.__types:
                    # Add unique station types (e.g. precip, runoff) if a 'Type' field exists in the metaheader
                    st = cnt  # last cnt becomes the starting column of the next type
                    order += 1

                # Information stored in __types array:
                # 1) Order that type was added in
                # 2) Starting index for data section
                # 3) Ending index for data section
                self.__types[words[self.__metaheader.index("Type")]] = [
                    order,
                    st,
                    cnt,
                ]
            except ValueError:
                if self.__verbose:
                    print('No "Type" metadata; skipping.')

            self.__stations.append(words[0])
            self.__stationIndex[words[0]] = cnt
            cnt += 1

        # Read the units and add to each type
        unittmp = next(it).split(":")[1].split(",")
        self.__headercount += 1
        for xx in unittmp:
            unit_pair = xx.split("=")
            self.__units[unit_pair[0].strip()] = unit_pair[1].strip()

        # Skip to the data section
        for line in it:
            self.__headercount += 1
            if line[0:10] == "##########":
                break

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Data section
        # The first 6 columns are [year month day hour minute seconds]
        thecols = ["year", "month", "day", "hour", "min", "sec"]

        # Add the remaining columns to the list
        for xx in self.__stations:
            thecols.append(xx)

        # Use pandas to read the data in from the remainder of the file
        # We use a custom date parser to convert the date information to a datetime
        # NOTE: 2023-03-21 skiprows option seems to be off by 1; test data starts
        #       at line 26, but skiprows=25 skips the first row of data.
        self.__rawdata = pd.read_csv(
            self.filename,
            skiprows=self.__headercount - 1,
            sep=r"\s+",
            header=0,
            names=thecols,
            engine="c",
            skipinitialspace=True,
            parse_dates={
                "time": ["year", "month", "day", "hour", "min", "sec"]
            },
            index_col="time",
        )

        self.__rawdata.index = pd.to_datetime(
            self.__rawdata.index, exact=True, cache=True, format=TS_FORMAT
        )

        # Convert the missing data (-999.0) to NaNs
        self.__rawdata.replace(
            to_replace=self.__missing, value=np.nan, inplace=True
        )

        self.__isloaded = True
