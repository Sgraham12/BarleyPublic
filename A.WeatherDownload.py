# Get Weather data
#Stations Sidney 2NW (1992-present), Lincoln 1500 N 45th (1989-present), Memphis 5N (Mead; 1994-present), Colby (Colby; 2010-present)

import requests
import pandas as pd
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
import csv
import json

class WebServices:
    def __init__(self, test=False):
        self.__server = 'https://awdn1.unl.edu/productdata/get?' if test else 'https://awdn2.unl.edu/productdata/get?'
        retry_strategy = Retry(total=5, status_forcelist=[429, 500, 502, 503, 504], backoff_factor=1)
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.__http = requests.Session()
        self.__http.mount("http://", adapter)
        self.__http.mount("https://", adapter)
        self.__http.headers['User-Agent'] = 'AWDN Web Services'

    def getList(self, product='scqc1440', network=None):
        params = {"list": product, "network": network} if network is not None else {"list": product}
        response = self.__http.get(self.__server, params=params)
        response.raise_for_status()
        try:
            return response.json()
        except ValueError as err:
            raise Exception(response.text)

    def getActive(self, product='scqc1440', network=None):
        params = {"active": product, "network": network} if network is not None else {"active": product}
        response = self.__http.get(self.__server, params=params)
        response.raise_for_status()
        try:
            return response.json()
        except ValueError as err:
            raise Exception(response.text)

    def getData(self, name, begin, end, product='scqc1440', format="csv", units='si', tz='CST', network="all", sensor="all"):
        params = {"name": name, "productid": product, "begin": begin, "end": end, "format": format}
        if units != 'si':
            params['units'] = units
        if tz != 'CST':
            params['tz'] = tz
        if network != 'all':
            params['network'] = network
        if sensor != 'all':
            params['sensor'] = sensor

        response = self.__http.get(self.__server, params=params)
        response.raise_for_status()
        if format == "csv":
            return response.text
        try:
            return response.json()
        except ValueError as err:
            raise Exception(response.text)

    def getGrid(self, date, product='scqc1440'):
        params = {"grid": product, "date": date}
        response = self.__http.get(self.__server, params=params)
        response.raise_for_status()
        try:
            return response.text
        except ValueError as err:
            raise Exception(response.text)

def flatten_row(row):
    """Flatten nested lists in a row dictionary."""
    flattened_row = {}
    for key, value in row.items():
        if isinstance(value, list):
            for i, item in enumerate(value):
                flattened_row[f"{key}_{i}"] = item
        else:
            flattened_row[key] = value
    return flattened_row

# Define the station and product parameters
stations = ["Sidney 2NW", "Lincoln 1500 N 45th", "Memphis 5N", "Colby"]  # List of stations
product_id = "scqc1440"

start_year = 2001
end_year = 2024

weather_service = WebServices(test=False)

for station_name in stations:
    all_data = []  # Reset data for each station
    for year in range(start_year, end_year + 1):
        begin_date = f"{year}0101"  # January 1st of the current year
        end_date = f"{year}1231"    # December 31st of the current year
        
        try:
            data_csv = weather_service.getData(name=station_name, begin=begin_date, end=end_date, product=product_id, format="csv")
            data = list(csv.DictReader(data_csv.splitlines()))
            all_data.extend(data)
            print(f"Data retrieved for {station_name} in {year}")
        except Exception as e:
            print(f"Error retrieving data for {station_name} in {year}: {e}")

    output_file = f"station_data_{station_name.lower().replace(' ', '_').replace(',', '').replace('/', '_')}.csv"

    if all_data:
        flattened_data = [flatten_row(row) for row in all_data]
        with open(output_file, "w", newline='') as csvfile:
            fieldnames = flattened_data[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(flattened_data)
        print(f"Data saved to {output_file}")
    else:
        print(f"No data retrieved for {station_name}.")

