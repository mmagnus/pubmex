#!/usr/bin/env python

# -*- coding: utf-8 -*-
import sys

import pygtk
pygtk.require("2.0")
import gtk
import gtk.glade
import pubmex 

class appGTK:
	def __init__(self):

		self.gladefile = "pubmex_gui.glade"
		self.wTree = gtk.glade.XML(self.gladefile) 

		self.window = self.wTree.get_widget("MainWindow")
		if (self.window):
			self.window.connect("destroy", gtk.main_quit)

		dic = { "on_get_button_clicked" : self.on_get_button_clicked}
		self.wTree.signal_autoconnect(dic)



	def on_get_button_clicked(self, widget):
		input=self.wTree.get_widget('input_text').get_text()#input_text
		status_bar=self.wTree.get_widget('status_bar')#fake status_bar

		if input:
			if pubmex.is_it_pmid(input):
				status_bar.set_label('connecting... PMID recognized')
				customed_title = self.wTree.get_widget('customed_title_entry').get_text().strip()
				output_text = pubmex.get_title_via_pmid(input, False, customed_title)
				output_buffer = self.wTree.get_widget('output_text').get_buffer()
				output_buffer.set_text(output_text)
				status_bar.set_label('done')
			else:
				status_bar.set_label('connecting... DOI recognized')
				customed_title = self.wTree.get_widget('customed_title_entry').get_text().strip()
				output_text = pubmex.get_title_via_doi(input, False, customed_title)
				output_buffer = self.wTree.get_widget('output_text').get_buffer()
				output_buffer.set_text(output_text)
				status_bar.set_label('done')
		else:
			status_bar.set_label('error: empty input, pass PMID or DOI!')

		

		
def set_status_bar(text):
	pass
		

if __name__ == "__main__":
	hwg = appGTK()
	gtk.main()
