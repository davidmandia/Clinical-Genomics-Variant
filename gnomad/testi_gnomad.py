import requests
from bs4 import BeautifulSoup

url = "https://gnomad.broadinstitute.org/variant/1-55524293-TG-T?dataset=gnomad_r2_1"

# Send a GET request to the URL
response = requests.get(url)

# Check if the response is successful
if response.status_code == 200:
    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.content, 'html.parser')
    
    # Extract specific data from the HTML content
    # Note: You'll need to inspect the webpage to determine the relevant tags and attributes to extract
    variant_data = {}
    
    # For example, extracting the title of the webpage
    title = soup.title.string
    variant_data['title'] = title
    
    # You can extract more specific data as needed by inspecting the HTML structure
    # For demonstration, extracting all text inside a specific tag (e.g., <h1>)
    # header = soup.find('h1').get_text() if soup.find('h1') else 'N/A'
    # variant_data['header'] = header
    
    print(soup)
else:
    print(f"Failed to retrieve data. HTTP Status Code: {response.status_code}")


