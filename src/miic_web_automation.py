#!/usr/bin/env python3
"""
MIIC Web Server Automation
Author: Ali Chemkhi

This module provides functions to automate the submission of datasets to the MIIC web server
at https://miic.curie.fr/workbench.php using Selenium WebDriver.
"""

import time
import os
from pathlib import Path
from typing import List, Optional, Tuple
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.chrome.options import Options


def submit_miic_jobs(
    csv_file_paths: List[str],
    st_txt_paths: List[str],
    job_names: Optional[List[str]] = None,
    email: str = "ali.chemkhi@curie.fr",
    headless: bool = False,
    wait_time: int = 10,
    orientation_cut: str = "0.01" , # Default value for orientation cut
    skeleton_cut: bool = True  # Default value for skeleton cut
) -> List[str]:
    """
    Submit multiple datasets to MIIC web server and return result links.
    
    Parameters:
    -----------
    csv_file_paths : List[str]
        List of paths to CSV files containing the datasets
    st_txt_paths : List[str]
        List of paths to corresponding .st.txt files
    job_names : Optional[List[str]]
        List of job names. If None, will use CSV file basenames
    email : str
        Email address for job notifications (default: ali.chemkhi@curie.fr)
    headless : bool
        Whether to run browser in headless mode (default: False)
    wait_time : int
        Time to wait between actions (default: 10 seconds)
    
    Returns:
    --------
    List[str]
        List of result URLs from MIIC web server
    """
    
    if len(csv_file_paths) != len(st_txt_paths):
        raise ValueError("Number of CSV files must match number of ST.txt files")
    
    # Generate job names if not provided
    if job_names is None:
        job_names = [Path(p).name.replace('.csv', '') for p in csv_file_paths]
    elif len(job_names) != len(csv_file_paths):
        raise ValueError("Number of job names must match number of CSV files")
    
    # Setup Selenium driver
    options = Options()
    options.headless = headless
    prefs = {"download.prompt_for_download": False}
    options.add_experimental_option("prefs", prefs)
    
    driver = webdriver.Chrome(
        service=Service(ChromeDriverManager().install()), 
        options=options
    )
    
    result_links = []
    
    try:
        for csv_path, st_path, name in zip(csv_file_paths, st_txt_paths, job_names):
            print(f"Processing: {name} - {csv_path}")
            
            try:
                # Step 1: Open MIIC workbench
                driver.get("https://miic.curie.fr/workbench.php")
                time.sleep(5)
                
                # Step 2: Fill job name
                name_input = driver.find_element(
                    By.XPATH, "/html/body/div[2]/form/div[1]/div/div[2]/div/input"
                )
                name_input.clear()
                name_input.send_keys(name)
                time.sleep(2)
                
                # Step 3: Fill email
                email_input = driver.find_element(
                    By.XPATH, "/html/body/div[2]/form/div[1]/div/div[3]/div/input"
                )
                email_input.clear()
                email_input.send_keys(email)
                time.sleep(2)
                
                # Step 4: Upload CSV file
                csv_upload = driver.find_element(
                    By.XPATH, "/html/body/div[2]/form/div[1]/div/div[4]/div[1]/input[2]"
                )
                csv_upload.send_keys(csv_path)
                time.sleep(2)
                
                # Step 5: Upload ST.txt file
                st_upload = driver.find_element(
                    By.XPATH, "/html/body/div[2]/form/div[2]/div[4]/div/div[3]/div[1]/input[1]"
                )
                st_upload.send_keys(st_path)
                time.sleep(2)

                # --- skeleton cut
                if skeleton_cut : 
                    # Wait until the radio button is present and clickable
                    radio_button = driver.find_element(By.XPATH, "/html/body/div[2]/form/div[2]/div[5]/div/div[1]/div/input[1]")
                    # Scroll it into view
                    driver.execute_script("arguments[0].scrollIntoView(true);", radio_button)
                    time.sleep(0.5)  # Give browser time to render
                    driver.execute_script("arguments[0].click();", radio_button)
                    # Confirm selection
                    print("Selected:", radio_button.is_selected())
        
                # --- orientation cut
                upload_input = driver.find_element(By.XPATH, "/html/body/div[2]/form/div[2]/div[6]/div/div/div/input")
                # Set the value directly using JS
                driver.execute_script("arguments[0].value = arguments[1];", upload_input, orientation_cut)
                time.sleep(1)



                
                # Step 6: Submit job
                print(f"Submitting MIIC job for {name}...")
                submit_button = driver.find_element(
                    By.XPATH, "/html/body/div[2]/form/div[2]/div[8]/div/div/input"
                )
                submit_button.click()
                time.sleep(wait_time)
                
                # Step 7: Get result link
                result_link_element = driver.find_element(
                    By.XPATH, "/html/body/div[2]/div/div[2]/div[3]/div[2]/div/table/tbody/tr[1]/td[2]/a"
                )
                result_link = result_link_element.get_attribute("href")
                result_links.append(result_link)
                print(f"----> Result link: {result_link}")
                
            except Exception as e:
                print(f"Error processing {name}: {str(e)}")
                result_links.append(f"ERROR: {str(e)}")
                
    finally:
        driver.quit()
    
    return result_links


def submit_single_miic_job(
    csv_file_path: str,
    st_txt_path: str,
    job_name: Optional[str] = None,
    email: str = "ali.chemkhi@curie.fr",
    headless: bool = False,
    wait_time: int = 10
) -> str:
    """
    Submit a single dataset to MIIC web server and return result link.
    
    Parameters:
    -----------
    csv_file_path : str
        Path to CSV file containing the dataset
    st_txt_path : str
        Path to corresponding .st.txt file
    job_name : Optional[str]
        Job name. If None, will use CSV file basename
    email : str
        Email address for job notifications (default: ali.chemkhi@curie.fr)
    headless : bool
        Whether to run browser in headless mode (default: False)
    wait_time : int
        Time to wait between actions (default: 10 seconds)
    
    Returns:
    --------
    str
        Result URL from MIIC web server
    """
    
    result_links = submit_miic_jobs(
        csv_file_paths=[csv_file_path],
        st_txt_paths=[st_txt_path],
        job_names=[job_name] if job_name else None,
        email=email,
        headless=headless,
        wait_time=wait_time
    )
    
    return result_links[0] if result_links else "ERROR: No result obtained"


if __name__ == "__main__":
    # Example usage
    csv_files = [
        "/Users/alichemkhi/Desktop/myProjects/2D_gast/output/v0.3/2D_gast_grouped_v0.3/2D_gast_grouped_v0.3_WNT_PS.376.csv",
        "/Users/alichemkhi/Desktop/myProjects/2D_gast/output/v0.3/2D_gast_grouped_v0.3/2D_gast_grouped_v0.3_WNT_aPS.383.csv"
    ]
    
    st_files = [
        "/Users/alichemkhi/Desktop/myProjects/2D_gast/output/v0.3/2D_gast_grouped_v0.3/2D_gast_grouped_v0.3_WNT_PS.376.st.txt",
        "/Users/alichemkhi/Desktop/myProjects/2D_gast/output/v0.3/2D_gast_grouped_v0.3/2D_gast_grouped_v0.3_WNT_aPS.383.st.txt"
    ]
    
    # Submit jobs
    results = submit_miic_jobs(csv_files, st_files)
    
    # Print results
    for i, result in enumerate(results):
        print(f"Job {i+1}: {result}")