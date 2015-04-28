__author__ = 'brunosousa'

import requests
import json
import numpy

class MeTH:
    def __init__(self, server_):
        self.server = server_
        self.url_post = '/ocpu/library/MeTH/R/run_METH_Json'

    def call_meth(self,  params):
        server_url_post = self.server + self.url_post

        hdr = {'Content-Type': 'application/json'}
        my_params = json.dumps(params)
        response = requests.post(server_url_post, headers=hdr, data=my_params)

        if (response.status_code == 201):
            aux_text = response.text.split('\n')[0] # get the val stored in .val

            server_url_get = self.server + aux_text + '/print'
            response = requests.get(server_url_get)

            return (response.json())
        else:
            return ''

    def call_meth_matrices(self,  mBen, mCost, wBen, wCost):
        server_url_post = self.server + self.url_post

        params={'mBen': str(mBen), 'mCost': str(mCost), 'weiBen': str(wBen), 'weiCost': str(wCost)}

        return(self.call_meth(params))


if __name__ == '__main__':
    server = 'http://localhost:7414'

    my_params = {'mBen': '[[1,2,3], [2,5,6], [3,1,1]]', 'mCost': '[[1,2,3], [2,5,6], [3,1,1]]', 'weiBen': '[0.5, 0.5]',
                  'weiCost': '[0.5,0.5]'}

    meth = MeTH(server)
    print(meth.call_meth(my_params))

    mBen = [[1,2,3],[2,5,6],[3,1,1]]
    mCost = [[1,2,3],[2,5,6],[3,1,1]]
    weiBen = [0.5, 0.5]
    weiCost = [0.5, 0.5]

    print(meth.call_meth_matrices(mBen,mCost,weiBen,weiCost))



