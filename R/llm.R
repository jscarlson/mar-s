#' @title Multi-Class Text Classifier using ellmer with OpenAI
#' @description Provides an interface to perform (mutually exclusive) multi-class text classification using the ellmer package
#' @import ellmer
#' @importFrom assertthat assert_that
#' @importFrom jsonlite toJSON

#' @param text Character vector containing the texts to classify
#' @param classes Character vector specifying the possible classes to classify against
#' @param examples Data frame with columns 'text' and 'class' containing example classifications.
#'                 Maximum of 10 rows allowed. Class values must be in the classes parameter.
#' @param model Character string specifying the model to use (default: "gpt-4.1")
#' @param max_tokens Integer specifying maximum tokens in response (default: 1000)
#' @param system_prompt Character string for system-level instructions (optional)
#' @param base_url Character string for the API endpoint (default: "https://api.openai.com/v1")
#' @param api_key Character string for the API key (default: uses OPENAI_API_KEY environment variable)
#' @param api_args Named list of additional arguments to pass to the API (optional)
#' @param echo One of "none", "output", or "all" to control console output (default: "none")
#' @param max_active Integer specifying maximum number of simultaneous requests (default: 10)
#' @param rpm Integer specifying maximum requests per minute (default: 3)
#' @return A data frame containing:
#'   \item{text}{The input text}
#'   \item{class}{The predicted class}
#'   \item{confidence}{Numeric value between 0 and 1 indicating the confidence in the classification}
#' @export
multiclass_text_classifier_openai <- function(
    text, classes, examples,
    model = "gpt-4.1", params = list(max_tokens = 100),
    system_prompt = NULL, base_url = "https://api.openai.com/v1",
    api_key = Sys.getenv("OPENAI_API_KEY"), api_args = list(),
    echo = "none", max_active = 1, rpm = 3) {
  
  # Input validation
  assertthat::assert_that(
    is.character(text),
    length(text) > 0,
    is.character(classes),
    length(classes) > 0,
    !any(duplicated(classes)),
    is.data.frame(examples),
    nrow(examples) > 0,
    nrow(examples) <= 10,
    all(c("text", "class") %in% names(examples)),
    all(examples$class %in% classes),
    is.character(model),
    is.list(params),
    is.null(system_prompt) || is.character(system_prompt),
    is.character(base_url),
    is.character(api_key),
    is.list(api_args),
    is.numeric(max_active),
    max_active > 0,
    is.numeric(rpm),
    rpm > 0
  )
  
  # Match echo argument
  echo <- match.arg(echo, choices = c("none", "output", "all"))
  
  # Define the classification type specification
  type_classification <- ellmer::type_object(
    "Multi-class classification with confidence score",
    class = ellmer::type_enum(
      "Predicted class",
      values = classes
    ),
    confidence = ellmer::type_number(
      "Confidence in the classification, ranging from 0.0 to 1.0"
    )
  )
  
  # Format examples into prompt
  example_prompt <- paste0("Here are some examples of multi-class classification:\n\n")
  for (i in seq_len(nrow(examples))) {
    example_prompt <- paste0(
      example_prompt,
      "Text: ", examples$text[i], "\n",
      "Class: ", examples$class[i], "\n\n"
    )
  }
  
  # Create prompts for each input text
  prompts <- lapply(text, function(t) {
    paste0(
      example_prompt,
      "Now classify this text into exactly one of the following classes. ",
      "Available classes: ", paste(classes, collapse = ", "), "\n",
      t
    )
  })
  
  # Initialize chat object with all parameters
  chat <- ellmer::chat_openai(
    system_prompt = system_prompt,
    base_url = base_url,
    api_key = api_key,
    model = model,
    params = ellmer::params(
      max_tokens = params$max_tokens
    ),
    api_args = api_args,
    echo = echo
  )

  print(prompts)
  
  # Make the API calls in parallel with structured output
  tryCatch({
    results <- ellmer::parallel_chat_structured(
      chat = chat,
      prompts = prompts,
      type = type_classification,
      max_active = max_active,
      rpm = rpm
    )
    
    # Convert results to data frame with input text
    data.frame(
      text = text,
      class = results$class,
      confidence = results$confidence,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    stop(sprintf("Error calling LLM: %s", e$message))
  })
}

#' @title Multi-Class Text Classifier using ellmer with Hugging Face
#' @description Provides an interface to perform (mutually exclusive) multi-class text classification using the ellmer package with Hugging Face models
#' @import ellmer
#' @importFrom assertthat assert_that
#' @importFrom jsonlite toJSON

#' @param text Character vector containing the texts to classify
#' @param classes Character vector specifying the possible classes to classify against
#' @param examples Data frame with columns 'text' and 'class' containing example classifications.
#'                 Maximum of 10 rows allowed. Class values must be in the classes parameter.
#' @param model Character string specifying the Hugging Face model to use (default: "meta-llama/Llama-3.1-8B-Instruct")
#' @param params Named list of model parameters (optional). Note that not all parameters may be supported by all models.
#' @param system_prompt Character string for system-level instructions (optional, note that some models may not support this)
#' @param api_key Character string for the API key (default: uses HUGGINGFACE_API_KEY environment variable)
#' @param api_args Named list of additional arguments to pass to the API (optional)
#' @param echo One of "none", "output", or "all" to control console output (default: "none")
#' @param max_active Integer specifying maximum number of simultaneous requests (default: 1)
#' @param rpm Integer specifying maximum requests per minute (default: 3)
#' @return A data frame containing:
#'   \item{text}{The input text}
#'   \item{class}{The predicted class}
#'   \item{confidence}{Numeric value between 0 and 1 indicating the confidence in the classification}
#' @export
multiclass_text_classifier_hf <- function(
    text, classes, examples,
    model = "meta-llama/Llama-3.1-8B-Instruct", params = NULL,
    system_prompt = NULL, api_key = Sys.getenv("HUGGINGFACE_API_KEY"),
    api_args = list(), echo = "none", max_active = 1, rpm = 3) {
  
  # Match echo argument
  echo <- match.arg(echo, choices = c("none", "output", "all"))
  
  # Input validation
  assertthat::assert_that(
    is.character(text),
    length(text) > 0,
    is.character(classes),
    length(classes) > 0,
    !any(duplicated(classes)),
    is.data.frame(examples),
    nrow(examples) > 0,
    nrow(examples) <= 10,
    all(c("text", "class") %in% names(examples)),
    all(examples$class %in% classes),
    is.character(model),
    is.null(system_prompt) || is.character(system_prompt),
    is.character(api_key),
    is.list(api_args),
    is.numeric(max_active),
    max_active > 0,
    is.numeric(rpm),
    rpm > 0
  )
  
  # Define the classification type specification
  type_classification <- ellmer::type_object(
    "Multi-class classification with confidence score",
    class = ellmer::type_enum(
      "Predicted class",
      values = classes
    ),
    confidence = ellmer::type_number(
      "Confidence in the classification, ranging from 0.0 to 1.0"
    )
  )
  
  # Format examples into prompt
  example_prompt <- paste0("Here are some examples of multi-class classification:\n\n")
  for (i in seq_len(nrow(examples))) {
    example_prompt <- paste0(
      example_prompt,
      "Text: ", examples$text[i], "\n",
      "Class: ", examples$class[i], "\n\n"
    )
  }
  
  # Create prompts for each input text
  prompts <- lapply(text, function(t) {
    paste0(
      example_prompt,
      "Now classify this text into exactly one of the following classes. ",
      "Available classes: ", paste(classes, collapse = ", "), "\n",
      "Please provide both the predicted class and a confidence score between 0 and 1.\n\n",
      "Text to classify: ", t
    )
  })
  
  # Initialize chat object with Hugging Face
  chat <- ellmer::chat_huggingface(
    system_prompt = system_prompt,
    model = model,
    params = ellmer::params(params),
    api_key = api_key,
    api_args = api_args,
    echo = echo
  )
  
  # Make the API calls in parallel with structured output
  tryCatch({
    results <- ellmer::parallel_chat_structured(
      chat = chat,
      prompts = prompts,
      type = type_classification,
      max_active = max_active,
      rpm = rpm
    )
    
    # Convert results to data frame with input text
    data.frame(
      text = text,
      class = results$class,
      confidence = results$confidence,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    if (grepl("API key", e$message, ignore.case = TRUE)) {
      stop("Invalid or missing API key. Please check your HUGGINGFACE_API_KEY environment variable.")
    } else if (grepl("rate limit", e$message, ignore.case = TRUE)) {
      stop("Rate limit exceeded. Try reducing the number of parallel requests or increasing the rpm parameter.")
    } else {
      stop(sprintf("Error calling Hugging Face API: %s", e$message))
    }
  })
}
