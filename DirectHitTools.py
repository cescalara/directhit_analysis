from ROOT import TFile
import numpy as np
from matplotlib import pyplot as plt

from ipywidgets import FloatProgress
from IPython.display import display

from contextlib import contextmanager
from collections import namedtuple

class DirectHitSearch():
    """
    Implement a direct hit search on EUSO data
    Input is the name of standard ROOT TFile
    """
    def __init__(self):

        # constants
        self._rows = 48
        self._cols = 48

        self.datafile = TFile() 
        self.n_gtu = 0
        self.pkt_len = 0

        # thresholds
        self.counts_threshold = 6
        self.duration_threshold = 2
        self.min_area = 10
        self.max_sum = 10e3
        
        # initialisation
        self.Events = namedtuple("Events",
                                "filename gtu time duration shape n_gtu n_event n_lines n_circles pkt_len")
        self.Events.filename = ""
        self.Events.gtu = []
        self.Events.time = []
        self.Events.duration = []
        self.Events.shape = []
        self.Events.n_gtu = 0
        self.Events.n_lines = 0
        self.Events.n_circles = 0
        self.Events.pkt_len = 128

        # settings
        self.set_progress = False
        self.set_analysis = True
        
    def print_search_params(self):
        """
        print the current search parameters
        """
        print 'Search paramters:'
        print '-----------------'
        print 'file: ' + self.datafile.GetName()
        print 'counts_threshold: ' + str(self.counts_threshold)
        print 'duration_threshold: ' + str(self.duration_threshold)
        print 'min_area: ' + str(self.min_area)
        print 'max_sum: ' + str(self.max_sum)
        print 'pkt_len: ' + str(self.pkt_len)
        
    @contextmanager
    def open(self, filename):
        """
        open a ROOT TFile for analysis
        * checks if the file have the tevent TTree
        * checks if the file is squeezed and sets self.pkt_len
        * uses context management to avoid leaked file descriptors
        """
        self.datafile = TFile(filename)
        
        # check if ROOT TTree written
        try:   
            self.n_gtu = self.datafile.tevent.GetEntries()
        except AttributeError:
            print 'No tevent TTree found in ', filename
            self.set_analysis = False
        
        # set Events filename 
        self.Events.filename = filename
        
        # check if squeezed file
        if "sqz-dis" in filename:
            self.pkt_len = 64
        else:
            self.pkt_len = 128

        # context     
        yield
        self.datafile.Close()
        

    def find_candidates(self):
        """
        Search through all gtu for events satisfying the following criteria:
        - signal is above threshold in at least 2 adjacent pixels

        returns an event list and a label list
        """
        from scipy import ndimage
        from itertools import groupby
        from operator import itemgetter

        pcd = np.zeros((1, 1, self._rows, self._cols), dtype = 'B')
        focal_surface = np.zeros((self._rows, self._cols), dtype = 'B')
        self.datafile.tevent.SetBranchAddress("photon_count_data", pcd)

        # display progress
        if self.set_progress == True:
            prog = FloatProgress(min = 0, max = self.n_gtu)
            display(prog)
        
        detection_gtu = []
        detection_label = []
        for k in range(self.n_gtu):
            self.datafile.tevent.GetEntry(k)
            focal_surface[:][:] = pcd[0][0][:][:]
       
            # find pixels above threshold
            blobs = focal_surface > self.counts_threshold 
            if blobs.any() == True:

                # label connected regions that satisfy this condition
                labels, nlabels = ndimage.label(blobs)
                # find labelled regions greater than a certain size
                size = np.bincount(labels.ravel())

                if np.max(size[1]) > self.min_area and np.sum(focal_surface) < self.max_sum:
                # store the event and its label
                    detection_gtu.append(k)
                    biggest_label = size[1:].argmax() + 1
                    detection_label.append(labels == biggest_label)

            if self.set_progress == True:
                prog.value += 1

        return detection_gtu, detection_label

    
    def rm_long_events(self, detection_gtu):
        """
        remove events longer than self.duration_threshold from an event list
        returns a filtered event list to the Events tuple
        """
        from itertools import groupby
        from operator import itemgetter

        # initialise
        self.Events.gtu = []
        self.Events.duration = []
        
        # find consecutive GTU runs in candidate events
        for key, group in groupby(enumerate(detection_gtu), lambda (i, x): i-x):
            gtu_range = map(itemgetter(1), group)
            len_gtu = len(gtu_range)

            if len_gtu <= self.duration_threshold:
                self.Events.gtu.append(gtu_range[0])
                self.Events.duration.append(len_gtu)

        return self.Events.gtu, self.Events.duration

    def classify_shape(self):
        """
        identify the shape of events in an event list
        return a list of the event shapes
        """
        from scipy import ndimage
        from skimage.measure import regionprops
        from itertools import groupby
        from operator import itemgetter

        # initialise
        pcd = np.zeros((1, 1, self._rows, self._cols), dtype = 'B')
        focal_surface = np.zeros((self._rows, self._cols), dtype = 'B')
        self.datafile.tevent.SetBranchAddress("photon_count_data", pcd)

        self.Events.shape = []
        for e in self.Events.gtu:

            self.datafile.tevent.GetEntry(e)
            focal_surface[:][:] = pcd[0][0][:][:]

            # find pixels above threshold
            blobs = focal_surface > self.counts_threshold 
            if blobs.any() == True:

                # label connected regions that satisfy this condition
                labels, nlabels = ndimage.label(blobs)

                # find the eccentricty of the largest object
                size = np.bincount(labels.ravel())
                biggest_label = size[1:].argmax() 
                props = regionprops(labels)
                ecc = props[biggest_label].eccentricity
                length = props[biggest_label].major_axis_length

                # if object is eccentric and long, classify as linear
                if ecc > 0.7 and length > 10:
                    self.Events.shape.append('linear')
                else:
                    self.Events.shape.append('circular')
                    
        return self.Events.shape
    
    def plot_focal_surface (self, gtu_num):
        """
        plot the focal surface for a given gtu_num
        """
        # initialise
        pcd = np.zeros((1, 1, self._rows, self._cols), dtype = 'B')
        self.datafile.tevent.SetBranchAddress("photon_count_data", pcd)
        focal_surface = np.zeros((self._rows, self._cols), dtype = 'B')

        # get the pixel values
        self.datafile.tevent.GetEntry(gtu_num)
        focal_surface[:][:] = pcd[0][0][:][:]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        p = plt.imshow(focal_surface, axes=ax, interpolation='nearest')
        plt.xlabel('pixel X')
        plt.ylabel('pixel Y')
        plt.title(self.datafile.GetName() + "\n" + 'GTU: ' + str(gtu_num), y = 1.1)

        # Set ticks
        major_ticks = np.arange(-.5, self._rows, 8)
        minor_ticks = np.arange(-.5, self._rows, 1)
        #ax.set_xlim(0.0,47.0)
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)
        ax.set_xticklabels(np.arange(0, self._rows + 1, 8));
        ax.set_yticklabels(np.arange(self._cols, 0, -8));

        # set grid
        #ax.grid(which='minor', alpha=0.2)
        ax.grid(color = 'k', linestyle = '-', linewidth = 2)
        ax.grid(which = 'major', alpha = 0.4)

        # add colourbar
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.8, 0.1, 0.05, 0.8])
        cbar = fig.colorbar(p, cax=cbar_ax)
        cbar.set_label('# of counts', labelpad = 1)
        cbar.formatter.set_powerlimits((0, 0))

        
    def add_file_summary(self):
        """
        create file summary and add to the Events tuple
        """
        self.Events.n_gtu = self.n_gtu
        self.Events.n_events = len(self.Events.gtu)
        self.Events.n_lines = self.Events.shape.count('linear')
        self.Events.n_circles = self.Events.shape.count('circular')
        self.Events.pkt_len = self.pkt_len

        
