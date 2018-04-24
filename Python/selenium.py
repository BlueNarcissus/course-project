from selenium import webdriver
from selenium.webdriver.common.keys import Keys

class journal(object):
    def __init__(self):
        self.driver=webdriver.Chrome()
        self.driver.get('http://www.nationalgeographic.com/')
    
    def print_content(self):
        content=self.driver.find_element_by_id('mostreadthelatest_lr_mostread')
        mostreads=self.driver.find_elements_by_class_name('mt3_fiveup-card-title')
        i=1
        for mostread in mostreads:
            print(str(i)+'. '+mostread.text+'\n')
            i=i+1
        self.quit() # quit browser completely!

    def quit(self):
        self.driver.quit()


new=journal()
new.print_content()
