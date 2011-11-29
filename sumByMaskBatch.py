import PyRaster
import os,sys,numpy,future
#--------------------------------------
def sumRegion(argList):
    mask=argList[0]
    data=argList[1]
    code=argList[2]
    
    data=numpy.where(mask==code,data,0.0)
    summa=numpy.sum(numpy.sum(data))
    return summa

maskFile="D:/SMED_temp/Geo_fordelning/data2005/Cat_02b/abcdeux2.asc"
wd="D:/SMED_temp/Geo_fordelning/data2005"
aggregateName="simairVedair"
units="gpers"
summaryFile="D:/SMED_temp/Geo_fordelning/data2005/Cat_02b/emis_tot_abcdeux.txt"
years=["2005"]#, "2010", "2020"]
substances=["NOx","PM10"]#"TA",
#substances=["NO2","SO2","PM10","NOx","CO","Bensen","NMVOC"]
sectors=["total"]#"residential","navigation","other"]#"traffic",
flippa=0
#--------------------
maskRast=PyRaster.raster()
maskRast.readAscii(maskFile)
maskRast.nodataToZero()
maskCodes=maskRast.unique()
if maskRast.nodata in maskCodes:
    maskCodes.remove(maskRast.nodata)
if 0.0 in maskCodes:
    maskCodes.remove(0.0)

rows=[]
for code in maskCodes:
    rows.append([code])
    
headerYear=[""]
headerSector=[""]
headerSubstance=["code"]
for sector in sectors:
    for year in years:
        for substance in substances:
            emisFile=os.path.join(wd,year,aggregateName,sector,sector+substance+units+year+".asc")
            print "EmisFile: ", emisFile
            emisRast=PyRaster.raster()
            
            try:
                emisRast.readAscii(emisFile)         
                print "Summing sector: ", sector, "year: ", year, "substance: ",substance
                
                headerYear.append(year)
                headerSubstance.append(substance)
                headerSector.append(sector)
                
                emisRast.nodataToZero()
                if flippa:
                    emisRast.flip()
                print "Finding all mask codes"
                 
                for i in range(len(maskCodes)):
                    ortRast=PyRaster.raster()
                    code=maskCodes[i]
                    #print "summing for mask code: ", code
                    sumThread=future.Future(sumRegion,[maskRast.data,emisRast.data,code])
                    sum=sumThread()
                    rows[i].append(sum)
            except:
                print "No ",substance, " in sector ", sector, " year: ",year

output=open(summaryFile,'w')            

for i in range(len(headerYear)):
    output.write(headerYear[i]+"\t")
output.write("\n")

for i in range(len(headerSector)):
    output.write(headerSector[i]+"\t")
output.write("\n")

for i in range(len(headerSubstance)):
    output.write(headerSubstance[i]+"\t")
output.write("\n")

for row in rows:
    for i in range(len(row)):
        output.write(str(row[i])+"\t")
    output.write("\n")

output.close()
print "finished!"
    