#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.


    To manage the update of the reference data, i.e. download, verification, and decompression.

    Directories are compressed using command: tar -cf - directory/ | xz -9 -c - > directory.tar.xz

    Files tar.xz are decompressed using command: tar -xJf directory.tar.xz
"""
import os
import urllib2
import hashlib
import subprocess
import time
import argparse

# block size when reading a file
BLOCK_SIZE = 8192

class Settings():
    def __init__(self):
        self._url = "http://algbio.cs.uni-duesseldorf.de/software/ppsp"
        self._customVal = ':'
        self._refVersion = None
        try:
            self._version = urllib2.urlopen(self._url + '/' + 'version.txt').read()
        except urllib2.HTTPError:
            self._version = "1_4"
            print("Can't get the current version from the server, version '%s' will be considered." % self._version)

    def getRemoteSrc(self):
        return self._url + "/" + self._version

    def getLocalDst(self, type):
        if type == 'ref':
            return "/mnt/host_shared"  # "/Users/ivan/Documents/nobackup/vm_ref_download_test"
        elif type == 'tools' or type == 'sys':
            return "/apps/pps"  # "/Users/ivan/Documents/nobackup/vm_ref_download_test"
        elif type == 'custom':
            return self._customVal.split(':')[1]
        else:
            raise Exception("Not supported destination type: %s" % type)

    def getFileName(self, type):
        if type == 'ref':
            return 'reference_' + self.getRefVersion() + '.tar.xz'
        elif type == 'tools':
            return 'tools.tar.xz'
        elif type == 'sys':
            return 'sys.tar.xz'
        elif type == 'custom':
            return self._customVal.split(':')[0]
        else:
            raise Exception("Not supported type: %s" % type)

    def getDataDescription(self):
        if type == 'ref':
            return 'Reference data'
        elif type == 'tools':
            return 'Tools data'
        elif type == 'sys':
            return 'System data'
        else:
            return 'Data'

    def getRefVersion(self):
        return self._refVersion

    def setCustomValues(self, customVal):
        self._customVal = customVal

    def setUrl(self, url):
        self._url = url

    def setVersion(self, version):
        self._version = version
        print("Version '%s' will be used!" % self._version)

    def setRefVersion(self, refVersion):
        self._refVersion = refVersion


def downloadFile(fileName, settings, type):
    """
        Downloads a file from a remote location and verifies it's checksum that is stored
        in a corresponding file named file.checksum.
    """
    urlPath = settings.getRemoteSrc() + '/' + fileName
    urlChecksum = urlPath + '.checksum'

    try:
        remoteFile = urllib2.urlopen(urlPath)
    except urllib2.HTTPError as e:
        print("File not found on the server: %s" % urlPath)
        raise e

    meta = remoteFile.info()
    fileSize = int(meta.getheaders("Content-Length")[0])

    print("Downloading file: %s " % urlPath)
    print("File size: %s MB" % round((int(fileSize) / (1024*1024)), 1))
    print("Progress:")

    localFilePath = os.path.join(settings.getLocalDst(type), fileName)
    try:
        localFile = open(localFilePath, 'wb')
    except Exception as e:
        print("\nCan't determine local destination for file: %s" % localFilePath)
        raise e

    fileRead = 0
    blockSize = BLOCK_SIZE
    try:
        while True:
            block = remoteFile.read(blockSize)
            if not block:
                break
            localFile.write(block)
            fileRead += len(block)
            percent = round((float(fileRead) / float(fileSize)) * 100.0, 1)
            print str(percent) + ' %\r',
        localFile.close()
    except Exception as e:
        print("An error occurred while a file has been downloaded.")
        raise e

    if fileRead == fileSize:
        print ("100.0 %")
        print("Download successful!")
    else:
        raise Exception("The downloaded file doesn't have the required size! (%s != %s)" % (fileRead, fileSize))

    print("Verifying the file checksum.")
    try:
        checksumFile = urllib2.urlopen(urlChecksum)
        checksumRemote = checksumFile.read()
    except Exception as e:
        print("Unable to determine remote file checksum.")
        raise e

    try:
        # checksumLocal = hashlib.md5(open(localFilePath).read()).hexdigest()
        checksumLocal = getChecksumForFile(localFilePath)
    except Exception as e:
        print("Unable to calculate the checksum of the local file")
        raise e

    if checksumLocal == checksumRemote.strip():
        print('Checksum OK!')
    else:
        raise Exception("Wrong checksum!")


def decompress(fileName, settings, type):
    """
        Decompresses a tar.xz file.
    """
    srcFilePath = os.path.join(settings.getLocalDst(type), fileName)

    if srcFilePath.split('.')[-1] != 'xz':
        raise("Trying to decompress a non-xz file: " + srcFilePath)

    cmd = str('tar -xJf ' + srcFilePath)
    print("Decompressing file: %s" % srcFilePath)
    try:
        cmdProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=os.path.dirname(srcFilePath), stdout=None)
        cmdProc.wait()
    except Exception as e:
        print("Can't run command: %s" % cmd)
        raise e

    if cmdProc.returncode != 0:
        print("Command returned with a non-zero status: %s: %s" % (cmdProc.returncode, cmd))
        raise Exception("Can't decompress file: %s" % srcFilePath)

    print("File '%s' was successfully decompressed." % fileName)


def getChecksumForFile(f):
    fr = open(f)
    md5 = hashlib.md5()
    while True:
        data = fr.read(BLOCK_SIZE)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()


def getChecksumForDir(dirPath):
    """
        For each file in the directory, enter the path relative to the dirname to a file and its checksum.
    """
    outFile = dirPath + '.checksum'
    out = open(outFile, 'w')
    checkSumWalk(dirPath, dirPath, out)
    out.close()

def checkSumWalk(f, root, out):
    if os.path.isfile(f):
        assert root in f
        if not os.path.basename(f).startswith('.'):
            out.write(os.path.sep.join(f[f.find(root.split(os.path.sep)[-1]):].split(os.path.sep)[1:]) +
                      '\t' + getChecksumForFile(f) + '\n')
            # + '\t' + hashlib.md5(open(f).read()).hexdigest() + '\n')
    else:
        for child in os.listdir(f):
            checkSumWalk(os.path.join(f, child), root, out)


def verifyChecksumForDir(dirName, settings, type):
    """
        For each file listed in a list, verifies its checksum.
    """

    dirName = os.path.join(settings.getLocalDst(type), dirName)

    urlChecksumFile = settings.getRemoteSrc() + '/' + os.path.basename(dirName) + '.checksum'
    print("Verification of the decompressed directory '%s' started." % dirName)
    try:
        for line in urllib2.urlopen(urlChecksumFile):
            f, checksum = line.split('\t')
            f = os.path.join(dirName, f)
            # if checksum.strip() != hashlib.md5(open(f).read()).hexdigest():
            if checksum.strip() != getChecksumForFile(f):
                raise Exception("File '%s' is corrupted, it has a wrong checksum." % f)
    except Exception as e:
        print("Unable to verify directory: %s" % dirName)
        raise e

    print("Checksum verification completed successfully!")


def updateData(settings, type):
    """
        Updates data in a folder, the existing folder is renamed.
    """
    print('%s preparation started!' % settings.getDataDescription())

    if not os.path.exists(settings.getLocalDst(type)):
        print("The destination directory '%s' doesn't exist." % settings.getLocalDst(type))
        return
    fileName = settings.getFileName(type)
    downloadFile(fileName, settings, type)

    dirName = os.path.join(settings.getLocalDst(type), fileName.rsplit('.', 2)[0])
    if os.path.isdir(dirName):
        newDirName = str(dirName + '_' + str(time.time()).split('.')[0])
        try:
            os.rename(dirName, newDirName)
            print("Directory '%s' already exists, it will be renamed to '%s'" % (dirName, newDirName))
        except Exception:
            print("Can't rename directory '%s', this directory will be overwritten." % dirName)
    decompress(fileName, settings, type)
    dirName = ".".join(fileName.split('.')[:-2])
    verifyChecksumForDir(dirName, settings, type)
    print("%s are ready to use!" % settings.getDataDescription())


def _main():
    parser = argparse.ArgumentParser(description='Updates data.', epilog='')

    parser.add_argument('-r', '--reference-data', nargs=1, help='Updates the reference data, enter which version', dest='r')

    parser.add_argument('-t', '--tools', action='store_true', help='Updates tools.', dest='t')

    parser.add_argument('-s', '--sys', action='store_true', help='Updates system data.', dest='s')

    parser.add_argument('-c', '--custom', nargs=1, help='file_name.tar.xz:local_dir_path', dest='c')

    parser.add_argument('-u', '--url', nargs=1, help='Use this url as a source of remote files.', dest='u')

    parser.add_argument('-v', '--version', nargs=1, help='Use this version.', dest='v')

    args = parser.parse_args()

    if not args.r and not args.t and not args.s and not args.c:
        parser.print_help()
        return

    settings = Settings()

    if args.u:
        settings.setUrl(args.u[0])

    if args.v:
        settings.setVersion(args.v[0])

    if args.s:
        updateData(settings, type='sys')

    if args.r:
        settings.setRefVersion(args.r[0])
        if not os.path.exists(settings.getLocalDst('ref')):
            print("The path to the shared folder is broken, verify that the shared folder was set correctly!")
        updateData(settings, type='ref')

    if args.t:
        updateData(settings, type='tools')

    if args.c:
        settings.setCustomValues(str(args.c[0]))
        updateData(settings, type='custom')


def _test():
    pass
    # getChecksumForDir('/Volumes/VerbatimSSD/work/vm_rel_1_3/reference_NCBI20121122')
    getChecksumForDir('/Volumes/VerbatimSSD/work/vm_1_4/tools')
    # getChecksumForDir('/Volumes/VerbatimSSD/work/vm_1_4/sys')
    # getChecksumForDir('/net/metagenomics/projects/PPSmg/release/1_4/nobackup/reference_NCBI20140513')

if __name__ == "__main__":
    _main()
    # _test()