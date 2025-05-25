## GenCoExplorer

GenCoExplorer is an API web to identify potential cell type markers and predict genes functions at cell type level in the context of Parkinson's Disease. It uses the database from [scCoExpNets](https://github.com/aliciagp/scCoExpNets).

## How to use it

To deploy the tool developed in this Final Degree Project, you must first access the GitHub repository where all the necessary data and code are stored.

Once in the repository, click on the **Code** button, which will display a dropdown menu with several options. You can either clone the repository or download a `.zip` file containing all the contents.

In any case, once the repository has been downloaded, navigate to the folder named `code`, In any case, once the repository has been downloaded, navigate to the folder named `QueryAPI.py` is located.

You must run this file, which will provide an address and a port that you can access through a web browser. For example, by copying the URL `http://127.0.0.1:5000`, into the browser's address bar, the tool’s main page will appear, and you will be able to start using it.

## First look

![Home page](PaginaPrincipal.png)
At the home page of the API, there're three types of queries searching by: gene symbol, GO term, or cell type. 
All of them have a similar structure, at the top there're headers for changing beetwen queries. At the left panel there's a `Run example` button that autocompletes the filters and the search bar. Here an example for a gene symbol query.

![Home page](Gene_Relevance_API2.png)
