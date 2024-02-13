#!/usr/bin/env python

import os
import time
from numpy import choose
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import pandas as pd
import csv


username= "your_username"
password= "your_password"
headless = True
timeout = 60
download_folder= "downloads"


def genesis_files():
    datesDF= pd.read_csv('datesDF-old.csv')
    # print(datesDF)
    for i in range(len(datesDF)):
        datesDF.loc[i,'processed'] = False
    with open('datesDF-meta.csv', 'w+') as f:
        datesDF.to_csv(f, index=False)
    datesDF= pd.read_csv('datesDF-meta.csv')
    download_folder= "downloads"

    gisaid_downloads= os.listdir(download_folder)
    with open('gisaid_database.fasta', 'w+') as f:
        f.write('1972-06-011972-06-30 header\n')


def logining_epicov(username, password, headless, timeout):

    print("logging in  epicov...")
    mimeType = "application/octet-stream,application/excel,application/vnd.ms-excel,application/pdf,application/x-pdf"
    profile = Options()
    service=  Service("./geckodriver")
    profile.set_preference("browser.download.folderList", 2)
    profile.set_preference("browser.download.manager.showWhenStarting", False)
    profile.set_preference("browser.download.dir", os.getcwd()+"/downloads")
    profile.set_preference("browser.helperApps.neverAsk.saveToDisk", mimeType)
    profile.set_preference("plugin.disable_full_page_plugin_for_types", mimeType)
    profile.set_preference("pdfjs.disabled", True)

    options = Options()
    driver = webdriver.Firefox(service=service, options=profile)

    driver.implicitly_wait(20)  
    wait = WebDriverWait(driver, timeout)

    driver.get('https://www.epicov.org/epi3/frontend#')
    waiting_sys_timer(wait)
    
    driver.find_element(By.NAME, 'login').send_keys(username)
    driver.find_element(By.NAME,'password').send_keys(password)
    driver.execute_script("return doLogin();")
    
   
    time.sleep(5)

    try:
        driver.execute_script('document.querySelector("iframe").remove()')
        driver.execute_script('document.querySelector(".sys_olcurtain").remove()')
        driver.execute_script('document.querySelector(".sys_curtain").remove()')
    except:
        print("iframe error")
    return driver

def find_date_intervals(driver, wait, initial_date, final_date):
    try:
        driver.execute_script("document.getElementById('sys_curtain').remove()")
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")

    except:
        print("find date intervals")
    driver.find_elements(By.CLASS_NAME,"sys-form-fi-date")[3]
    time.sleep(5)
    try:
        driver.execute_script("document.getElementById('sys_curtain').remove()")
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")

    except:
        print("clicked on the from date box")
    #d river.switch_to.default_content()
    #WebDriverWait(driver, 30).until(EC.visibility_of_element_located((By.CLASS_NAME, 'ui-datepicker-calendar')))

    #EC.visibility_of_element_located(By.CLASS_NAME, 'ui-datepicker-calendar')
    # datefield = driver.find_elements(By.CLASS_NAME,"sys-form-fi-date")[0]
    time.sleep(wait)
    try:
        driver.execute_script("document.getElementById('sys_curtain').remove()")
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")

    except:
        print("line exception, didnt find the  sysolcurtain")
    # datefield = driver.find_elements(By.CLASS_NAME,"sys-form-fi-date")[3]
    time.sleep(5)
    ActionChains(driver).move_to_element(driver.find_elements(By.CLASS_NAME,"sys-form-fi-date")[0]).click().send_keys(Keys.COMMAND + "a").send_keys(Keys.DELETE).perform()

    ActionChains(driver).move_to_element(driver.find_elements(By.CLASS_NAME,"sys-form-fi-date")[0]).click().send_keys(initial_date).perform()
    time.sleep(5)
    # ActionChains(driver).move_to_element(driver.find_elements(By.CLASS_NAME,"sys-form-fi-date")[0]).click().send_keys(initial_date).perform()

    # driver.find_element(By.XPATH,'//*[@id="ce_rtopib_ak_input"]').clear()
    # driver.find_element(By.XPATH,'//*[@id="ce_rtopib_ak_input"]').send_keys(initial_date)
    time.sleep(wait)
    try:
        driver.execute_script("document.getElementById('sys_curtain').remove()")
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")

    except:
        print("set initial date")
    driver.find_element(By.XPATH,'//*[@id="ce_rts4fd_aq_input"]').click()
    try:
        driver.execute_script("document.getElementById('sys_curtain').remove()")
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
    except:
        print("settle input")
    time.sleep(10)
    datefield = driver.find_elements(By.CLASS_NAME,"sys-form-fi-date")[1]

    ActionChains(driver).move_to_element(datefield).click().send_keys(final_date).perform()
    driver.find_element(By.XPATH,'//*[@id="ce_rts4fd_aq_input"]').clear()
    driver.find_element(By.XPATH,'//*[@id="ce_rts4fd_aq_input"]').send_keys(final_date)
    
    time.sleep(wait+2)
    driver.execute_script("document.getElementById('sys_timer').remove()")


    driver.find_element(By.XPATH,'/html/body/form/div[3]/div/div[2]/div/div[3]/div[2]/div/div[2]/div/button').click()
    # search_btn = driver.find_element(By.XPATH,'/html/body/form/div[5]/div/div[2]/div/div[3]/div[4]/div/div[3]/div/button')
    # ActionChains(driver).move_to_element(search_btn).click().click().perform()
    time.sleep(15)
    #WebDriverWait(driver, 30).until(EC.visibility_of_element_located((By.CSS_SELECTOR, '#c_rcvys8_jj_leftinfo > span:nth-child(1)')))
    resultField = driver.find_element(By.CLASS_NAME,"sys-datatable-info-left")

    print(resultField.text)
    if resultField.text == "Total: 0 isolates":
        driver.find_element(By.CSS_SELECTOR,'#ce_rcvys8_ds > div:nth-child(1) > button:nth-child(1)').click()
        time.sleep(2)
        return False
    else :
        return True

def downloadFastaAndCloseDriver(driver):
    try:
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
        driver.execute_script("document.getElementById('sys_curtain').remove()")

    except:
        print("line 113")
    
    driver.find_elements(By.XPATH, "//input[@type='checkbox']")[0].click()
    time.sleep(20)
    try:
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
        driver.execute_script("document.getElementById('sys_curtain').remove()")

    except:
        print("line 120")
    time.sleep(20)
    resultField = driver.find_element(By.XPATH,'//div[@class="sys-datatable-info-left"]/span[1]')
    print(resultField.text)
    time.sleep(30)
    if resultField.text == "Total: 0 isolates":
        return True
    try:
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
        driver.execute_script("document.getElementById('sys_curtain').remove()")

    except:
        print("line 115")
        
    driver.find_elements(By.XPATH,'/html/body/form/div[5]/div/div[2]/div/div[3]/div[2]/div/div[3]/div/button')[0].click()
    # time.sleep(2)               /html/body/form/div[5]/div/div[2]/div/div[3]/div[2]/div/div[3]/div/button
    # driver.find_element(By.XPATH,'//button[@class="sys-event-hook sys-form-button"]')[2].click()
    time.sleep(30)
    try:
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
        driver.execute_script("document.getElementById('sys_curtain').remove()")

    except:
        print("line 125")

    iframe = driver.find_element(By.XPATH,"//iframe")
    driver.switch_to.frame(iframe) 





    # time.sleep(5)
    # driver.find_elements(By.XPATH, "//input[@value='all']")[0].click()
    # time.sleep(5)
    # try:
        
    #     driver.execute_script("document.getElementById('sys_olcurtain').remove()")
    #     driver.execute_script("document.getElementById('sys_curtain').remove()")
        
    # except:
    #     print("line 150")
    #     time.sleep(10)
    driver.find_elements(By.XPATH,'//button[@class="sys-event-hook sys-form-button"]')[2].click()
    time.sleep(30)
    ##AND GO BACK
    try:
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
    except:
        print("line 151")
    # time.sleep(20)
    # driver.find_elements(By.XPATH,'//button[@class="sys-event-hook sys-form-button"]')[0].click()
    # time.sleep(10)
    ##AND GO BACK
    driver.switch_to.default_content()
    time.sleep(30)
    try:
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
    except:
        print("line 159")
    time.sleep(15)
    WebDriverWait(driver, 400).until(EC.invisibility_of_element_located((By.CLASS_NAME,'sys_olcurtain')))

    print("downloading...")
    download_status=True
    while download_status:
        part_found = False
        for i in os.listdir("downloads"):
            if i.endswith(".part"):
                print("...")
                time.sleep(10)
                part_found = True
                break
        if not part_found:
            print("done")
            download_status=False
            break
    time.sleep(30)
    try:
        WebDriverWait(driver, 250).until(EC.invisibility_of_element_located((By.CLASS_NAME,'sys_olcurtain')))
    except:
        print("driver refresh")
        driver.refresh()
        WebDriverWait(driver, 10).until(EC.alert_is_present())
        driver.switch_to.alert.accept()  
        time.sleep(5)   
        
    WebDriverWait(driver, 250).until(EC.invisibility_of_element_located((By.CLASS_NAME,'sys_olcurtain')))
    # driver.find_elements(By.XPATH,'//button[@class="sys-event-hook sys-form-button"]')[0].click()
    driver.quit()
    return True
  
def download_fasta(driver, wait, accession_id):

    print(f'downloading {accession_id}...')
    input = driver.find_element(By.XPATH,"//div[@class='sys-form-fi-entry']//input")

    input.send_keys(accession_id)
    waiting_sys_timer(wait)

    tr = driver.find_element(By.XPATH,"//tbody/tr[@class='yui-dt-rec yui-dt-first yui-dt-last yui-dt-even']")
    ActionChains(driver).double_click(tr).perform()
    waiting_sys_timer(wait, 15)
    
    driver.execute_script("window.scrollTo(0,document.body.scrollHeight)")
    iframe = driver.find_element(By.XPATH,"//iframe")
    driver.switch_to.frame(iframe) 
    
    driver.find_element(By.XPATH,"/html/body/form/div[5]/div/div[2]/div[2]/div/div[2]/div/button").click()
    waiting_sys_timer(wait, 5)
    driver.find_element(By.XPATH,"/html/body/form/div[5]/div/div[2]/div[2]/div/div[1]/div/button").click()
    
    driver.switch_to.default_content()
    input.clear()
    waiting_sys_timer(wait, 5)

    return driver, wait

def waiting_sys_timer(wait, sec=1):
    wait.until(EC.invisibility_of_element_located((By.XPATH,  "//div[@id='sys_timer']")))
    time.sleep(sec)

def search_download_and_close_driver(driver, wait, initial_date, final_date):

    time.sleep(5)
    responseFromIntervalSearch=find_date_intervals(driver, 2, initial_date, final_date)
    time.sleep(10)

    if responseFromIntervalSearch:
        try:
            driver.execute_script("document.getElementById('sys_olcurtain').remove()")
            driver.execute_script("document.getElementById('sys_curtain').remove()")

        except:
            print("line 204")
        resultField = driver.find_element(By.XPATH,'//div[@class="sys-datatable-info-left"]/span[1]')
        print(resultField.text)
        if resultField.text != "Total: 0 isolates":
            print("line 300")
            downloadFastaAndCloseDriver(driver)
            time.sleep(10)
        
    else : 
        print("No data found")
        return False
    
def row_to_process(datesDF):
    for i in range(len(datesDF)):
        if datesDF.loc[i,'processed'] == False:
            print(datesDF.loc[i,'initial_date'], datesDF.loc[i,'final_date'])
            if (datesDF.loc[i,'initial_date']=="end"):
                print ("all dates have been processed")
                exit()
            return i
    

    print("All dates have been processed")
    return -1
    
# driver= logining_epicov(username, password, headless, timeout)
# driver.find_element(By.LINK_TEXT,"EpiFlu™").click()
# time.sleep(10)

# driver.execute_script("return sys.call('c_rtopib_h9','GoSearch',{});")
# # driver.find_element(By.XPATH,"/html/body/form/div[5]/div/div[2]/div/div[2]/div/div/div/div[7]/div").click()
# time.sleep(20)

#load the data frame information
datesDF= pd.read_csv('datesDF-meta.csv')
#save download files state
res =row_to_process(datesDF) #find the next date to process



#if the response is -1 then the dates have been processed
while res != -1:
#set the arguments for the processing dates
    datesDF= pd.read_csv('datesDF-meta.csv')
    res =row_to_process(datesDF) #find the next date to process until it returns -1
    initial_files_at_download= os.listdir(download_folder)
    initial_date= datesDF.loc[res,'initial_date']
    final_date= datesDF.loc[res,'final_date']
    driver= logining_epicov(username, password, headless, timeout)
    driver.find_element(By.LINK_TEXT,"EpiFlu™").click()
    time.sleep(10)

    driver.execute_script("sys.call('c_rts4fd_hf','GoSearch',{});")
    # driver.find_element(By.XPATH,"/html/body/form/div[5]/div/div[2]/div/div[2]/div/div/div/div[7]/div").click()
    time.sleep(20)

    try:
        driver.execute_script("document.getElementById('sys_curtain').remove()")
        driver.execute_script("document.getElementById('sys_olcurtain').remove()")
    except:
        print("getting initial dates and listing files to download")
    try: 
        resultField = driver.find_element(By.XPATH,'//div[@class="sys-datatable-info-left"]')
    except: 
        resultField = driver.find_element(By.XPATH,'//div[@class="sys-form-fi-info"]')

    print(resultField.text)
    if resultField.text == "Total: 0 isolates":
        driver.find_element(By.XPATH, "/html/body/form/div[4]/div/div[2]/div/div[3]/div[2]/div/div[1]/div/button").click()
        datesDF= pd.read_csv('datesDF-meta.csv')
        res =row_to_process(datesDF) #find the next date to process
        initial_files_at_download= os.listdir(download_folder)
        initial_date= datesDF.loc[res,'initial_date']
        final_date= datesDF.loc[res,'final_date']
        try:
            driver.execute_script("document.getElementById('sys_curtain').remove()")
            driver.execute_script("document.getElementById('sys_olcurtain').remove()")
        except:
            print("line 226")
        # driver.find_element(By.XPATH,'/html/body/form/div[4]/div/div[2]/div/div[3]/div[2]/div/div[5]/div[1]/div').click()
        # time.sleep(2)
    time.sleep(25)

    data_found= search_download_and_close_driver(driver, 2, initial_date, final_date)


    time.sleep(15)
    if data_found:
    #find the newly downloaded file and change its name
        final_files_after_download= os.listdir(download_folder)
        difference2 = [element for element in final_files_after_download if element not in initial_files_at_download]
        if len(difference2)>0:
            os.rename(download_folder+"/"+difference2[0], download_folder+"/"+initial_date+final_date+".fasta")
            #if this file can successfully change its name then:
            #change the state of the date in the data frame and record it as a backup
            datesDF.loc[res,'processed']= True
            with open('datesDF-meta.csv', 'w+') as f:
                datesDF.to_csv(f, index=False)
    else:
        datesDF.loc[res,'processed']= True
        with open('datesDF-meta.csv', 'w+') as f:
            datesDF.to_csv(f, index=False)
        

exit()

