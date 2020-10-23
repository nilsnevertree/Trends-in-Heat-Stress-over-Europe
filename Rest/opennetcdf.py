#%%
# import sys
# import time
# import urllib.request
# import numpy as np
# import xarray as xr
# import os

# def get_paths(originpath = False, savepath = False):
#     if not originpath :
#         done = False
#         while not done: 
#             originpath = input('type in originpath\n')
#             print(originpath)

#             if input('correct filepath? y/n\n') == 'y':
#                 done = True
#     if not savepath :
#         path = os.getcwd()
#         if input('correct folderpath? \n{}\n y/n\n'.format(path)) == 'y':
#             done = True
#         else :
#             done = False
#         while not done: 
#             path = input('type in filepath\n')
#             print(path)
#             if input('correct path? y/n\n') == 'y':
#                 done = True

#         done = False
#         while not done: 
#             savename = input('type in savename\n')
#             savepath = path  + '\\' + savename
#             print(savepath)
#             if input('correct savename? y/n\n') == 'y':
#                 done = True

#     return originpath, savepath

# def reporthook(count, block_size, total_size):
#     global start_time
#     if count == 0:
#         start_time = time.time()
#         return
#     duration = time.time() - start_time
#     progress_size = int(count * block_size)
#     speed = int(progress_size / (1024 * duration))
#     percent = count * block_size * 100 / total_size
#     sys.stdout.write("\r...{:.2f}%, {:.1f} MB, {} KB/s, {} seconds passed".format(percent, progress_size / (1024 * 1024), speed, int(duration)) )
#     sys.stdout.flush()

# def save(url = False, filename = False):
#     url, filename = get_paths(originpath = url, savepath = filename)
#     urllib.request.urlretrieve(url, filename, reporthook)

# save()

# print('==========\ndone')
# input('press any key')

#%%
import requests
import sys
import time
url = "https://ftp.geomar.de/users/swahl/FOCI1.22-SW132_echam6_tracer_1850-2013_O3_pl_monthly_99615-1_zonmean.nc"

#Get the headers of the remote file
h = requests.head(url, allow_redirects=True)

#Get the size of the file
total_size = int(h.headers.get('content-length'))

#Request the file download with stream set to True
r = requests.get(url, stream=True)

#Open a local file for writing
localfile = open("mycopy2.nc", "wb")
chunks = 0

global global_start_time
global_start_time = time.time()
start_time = time.time()
#Process the file as it arrives
for chunk in r.iter_content(chunk_size=512):
    if chunk:

        chunks += 1
        downloaded = chunks * 512
        # An approximation as the chunks don't have to be 512 bytes
        progress = (downloaded/total_size)*100
        # print ("Download Progress",str(int(progress)),"%")
        duration = start_time - time.time()
        speed = 512
        # if duration != 0:
        #     speed = 512/duration
        # else:
        #     speed = 0
        sys.stdout.write("\r...{:.2f}%, {:.1f} MB, {:.1f} KB/s, {:.1f} seconds passed".format(progress, downloaded / (1024 * 1024), speed, time.time() - global_start_time))
        sys.stdout.flush()
        localfile.write(chunk)
    start_time = time.time()
print("Finished")
import requests
import sys
import time

def downloadFile(url, directory) :
  localFilename = url.split('/')[-1]
  with open(directory + '/' + localFilename, 'wb') as f:
    start = time.clock()
    r = requests.get(url, stream=True)
    total_length = r.headers.get('content-length')
    dl = 0
    if total_length is None: # no content length header
      f.write(r.content)
    else:
      for chunk in r.iter_content(1024):
        dl += len(chunk)
        f.write(chunk)
        done = int(50 * dl / total_length)
        sys.stdout.write("\r[%s%s] %s bps" % ('=' * done, ' ' * (50-done), dl//(time.clock() - start)))
        print ''
  return (time.clock() - start)

def main() :
  if len(sys.argv) > 1 :
        url = sys.argv[1]
  else :
        url = raw_input("Enter the URL : ")
  directory = raw_input("Where would you want to save the file ?")

  time_elapsed = downloadFile(url, directory)
  print "Download complete..."
  print "Time Elapsed: " + time_elapsed


if __name__ == "__main__" :
  main()