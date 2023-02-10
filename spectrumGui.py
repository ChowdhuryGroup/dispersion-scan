# -*- coding: utf-8 -*-

"""

GUI for Spectrometer

UPDATED: 12-04-2019



 M easurement

 A cquisition

'N

 D ata

 O rganization to

 R eplace

 L abVIEW

 A pplications



"""

import oceanOpticSpectrosco as s

import numpy as np
import tkinter as tk


import matplotlib

#matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from matplotlib.figure import Figure

import matplotlib.animation as ani

from matplotlib import style

style.use('ggplot')



def end():

    m.quit()

    m.destroy()

    return



def setmono():

    mono.comset(COMport.get())

    minfo.set(mono.info())

    monostate.set(mono.state())

    return 



def getcomnum():

    COMport.set(comnum.get())

    print(COMport.get())

    setmono()

    

    return



def setwave():

    wl.set(wlnum.get())

    mono.setwl(wl.get())

    monostate.set(mono.state())

    return



def setgrat():

    gr.set(gratnum.get())

    mono.setgr(gr.get())

    monostate.set(mono.state())

    return







def listdev():

    ools=s.sb.list_devices()

    lst=str(ools).replace(',','\n')

    lst=lst[1:-1:]

    oolist.set(lst)

    return



def averagespec(sernum,num):

    i=1

    avgspec=s.ocean(sernum).getspec()/num

    while i<num:

        avgspec=avgspec+s.ocean(sernum).getspec()/num

        i=i+1

    return avgspec



def savedata():

    i=1

    num=int(filnum.get())

    filename=tk.filedialog.asksaveasfilename(confirmoverwrite=True,title='MANDORLA - Choose Folder and Base File Name')

    while i<=num:

        visdata=OOsp(visspec).save()

        swirdata=OOsp(swirspec).save()

        filename=filename.replace('.txt','')

        

        if i<10:

            index='0'+str(i)

        else:

            index=str(i)

            

        np.savetxt(filename+'_vis_'+index,visdata, delimiter='\t', header='')

        np.savetxt(filename+'_swir_'+index,swirdata, delimiter='\t', header='')

        print('Saving file index: '+index)

        

        a=ooT1.get()

        b=ooT2.get()

        if a<b:

            pause=b

        else:

            pause=a

        

        s.time.sleep(pause/1000) # sleep is in seconds

        i=i+1

        

    print ('Save Complete!')    

    return



def OOselect():

    f1.setsernum(oonum.get())

    return



def animate(i,fig):

     

        plotsp=s.ocean(fig.sernum).getspec()

        

        x=plotsp[0,...]

        y=plotsp[1,...]

        

        

        fig.g.clear()

        fig.graph()

        fig.g.plot(x,y)

        return plotsp



def pullparam():

    

    return



def settimes():

    num1=ooT1.get()

    num2=ooT2.get()

    num3=ooT3.get()

    

    OOsp(visspec).setint(num1)

    OOsp(swirspec).setint(num2)

    OOsp('USB4F11156').setint(num3)

    return



class OOsp:

    

    def __init__(self, sernum):

        self.sernum=sernum

        self.f=Figure(figsize=(6,4),dpi=100)

        self.g=self.f.add_subplot(111)

        return

    

    def graph(self):

        self.g.set_title(self.sernum)

        self.g.set_xlabel('Wavelength (nm)')

        self.g.set_ylabel('Counts')

        return 



    def save(self):

        data=s.ocean(self.sernum).getspec()

        return data

    

    def setsernum(self,num):

        self.sernum=num

        return

    

    def setint(self,num):

        s.ocean(self.sernum).setinttime(num)

        return

 

visspec='OFX01123'

swirspec='NQ5200084'    

    

#make GUI window

m=tk.Tk()

m.title('MANDORLA')

rowNum=10;

colNum=3;



#frames

filefr=tk.Frame(m)

filefr.grid(row=1,column=0)

monofr=tk.Frame(m)

monofr.grid(row=1,column=1)

vistoolfr=tk.Frame(m)

vistoolfr.grid(row=7,column=1)

oofr=tk.Frame(m)

oofr.grid(row=8,column=0)





mono=s.mono()

COMport=tk.IntVar()

minfo=tk.StringVar()

minfo.set('    ****NOT CONNECTED****    ')

wl=tk.DoubleVar()

gr=tk.IntVar()

monostate=tk.StringVar()

monostate.set('Wavelength:  ------  || Grating: --  ---- g/mm BLZ=  ---NM')

oolist=tk.StringVar()

listdev()

oonum=tk.StringVar()

oonum.set(visspec)

ooT1=tk.IntVar()

ooT1.set(100)

ooT2=tk.IntVar()

ooT2.set(100)

ooT3=tk.IntVar()

ooT3.set(100)



f1=OOsp(oonum.get())



Visplot=FigureCanvasTkAgg(f1.f,m)

Visplot.draw()

Visplot.get_tk_widget().grid(row=8,column=1)

Vistool=NavigationToolbar2Tk(Visplot,vistoolfr)

Vistool.update()









# Save file widget

tk.Label(filefr,text='Number of Files to Save:').grid(row=0,column=0)

filnum=tk.Entry(filefr,width=5)

filnum.grid(row=0,column=1,columnspan=2)

saveb=tk.Button(filefr, text='Save', command=savedata)

saveb.grid(row=0,column=3)



#set COM port line

tk.Label(monofr,text='COM port:').grid(row=1,column=0)

comnum=tk.Spinbox(monofr,from_=1,to=9,width=10)

comnum.grid(row=1,column=1)

setcom=tk.Button(monofr,text='Set COM port #',command=getcomnum)

setcom.grid(row=1,column=2)                 



#set wavelength line

tk.Label(monofr,text='Wavelength (nm):').grid(row=3,column=0)                 

wlnum=tk.Entry(monofr, width=5)

wlnum.grid(row=3,column=1)

setWL=tk.Button(monofr,text='Set Wavelength',command=setwave)

setWL.grid(row=3,column=2)



#set grating line

tk.Label(monofr,text='Grating # :').grid(row=4,column=0)

gratnum=tk.Spinbox(monofr,from_=1,to=3, width=5)

gratnum.grid(row=4,column=1)

setg=tk.Button(monofr,text='Set Grating',command=setgrat).grid(row=4,column=2)

         

#status line

monostatus=tk.Label(monofr,textvariable=monostate)

monostatus.grid(row=5,column=0, columnspan=3)

monostatus.config(relief='ridge',pady=5,padx=5)





#make CLOSE button

stop=tk.Button(m,text='Exit', width=25, command=end)

stop.grid(row=rowNum,column=1)



#serial and model info

tk.Label(monofr,text='Monochromator:').grid(row=2,column=0)

minfobox=tk.Label(monofr,textvariable=minfo)

minfobox.grid(row=2,column=1, columnspan=2)

minfobox.config(relief='ridge', pady=5,padx=5)



#list OceanOptics spectrometers

tk.Label(oofr,text='OceanOptics specs:').grid(row=5, columnspan=2)

updateOO=tk.Button(oofr,text='Update List', command=listdev)

updateOO.grid(row=4, columnspan=2)

OOlst=tk.Label(oofr,textvariable=oolist)

OOlst.grid(row=6, columnspan=2)



tk.Label(oofr,text='Plot:').grid(row=0,column=0)

oosel1=tk.Radiobutton(oofr, text= 'FlameX', variable=oonum,value=visspec, command=OOselect)

oosel1.grid(row=1,column=0)

oosel2=tk.Radiobutton(oofr, text= 'NIRQuest', variable=oonum,value=swirspec, command=OOselect)

oosel2.grid(row=2,column=0)

oosel3=tk.Radiobutton(oofr, text= 'USB4000', variable=oonum,value='USB4F11156', command=OOselect)

oosel3.grid(row=3,column=0)



setT=tk.Button(oofr,text='Set Int. Time (ms)',command=settimes)

setT.grid(row=0,column=1,columnspan=2)

ootime1=tk.Entry(oofr,width=5, textvariable=ooT1)

ootime1.grid(row=1,column=1)

ootime2=tk.Entry(oofr,width=5,textvariable=ooT2)

ootime2.grid(row=2,column=1)

ootime3=tk.Entry(oofr,width=5, textvariable=ooT3)

ootime3.grid(row=3,column=1)





visp=ani.FuncAnimation(f1.f,animate,fargs=(f1,),interval=200)





m.mainloop()



