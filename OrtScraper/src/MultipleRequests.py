# -*- coding: utf-8 -*-
"""

"""

import grequests
import requests

class MultipleRequests:
    def __init__(self, urls, size):
        self.urls = urls
        self.size = size

    def exception(self, request, exception):
        print("Problem: {}: {}".format(request.url, exception))

    def async(self):
        session = requests.Session()
        results = grequests.map((grequests.get(u, stream=False, session=session) for u in self.urls),
                                exception_handler=self.exception, size=self.size, gtimeout=60)
        return results
