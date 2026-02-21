library(torch)
#library(luz)
#library(torchvision)
#library(torchdatasets)

# dataset
testDataset = dataset("test",
                      initialize = function(data)
                      {
                        self$x = torch_tensor(data$x)
                        self$y = torch_tensor(data$y)
                      },
                      .getitem = function(i) 
                      {
                        return(list(x = self$x[i, drop = FALSE], y = self$y[i, drop = FALSE]))
                      },
                      
                      .length = function() {
                        return(self$y$size()[[1]])
                      })

testData = tibble(x = seq(1, 50, by = 0.1), y = 1.37 + 2.4201*x^0.6516)
train_indices <- sample(1:nrow(testData), 0.5 * nrow(testData))

train_ds <- testDataset(testData[train_indices, ])
valid_ds <- testDataset(testData[setdiff(1:nrow(testData), train_indices), ])

train_dl <- train_ds %>% dataloader(batch_size = 8, shuffle = TRUE)
valid_dl <- valid_ds %>% dataloader(batch_size = 8, shuffle = FALSE)

net <- nn_module("linear",
                 initialize = function() 
                 {
                   self$m <- nn_linear(in_features = 1, out_features = 1)
                 },
                 forward = function(x) 
                 {
                   x %>% self$m()
                 })

trainingStart = Sys.time()
fitted <- net %>% setup(loss = function(y_hat, y_true) { return(nnf_mse_loss(y_hat, y_true)) }, # can this be simplified to just loss = nnf_mse_loss?
                        optimizer = optim_adam) %>%
  set_opt_hparams(lr = 0.1) %>%
  fit(train_dl, epochs = 100, valid_data = valid_dl)
Sys.time() - trainingStart

fitted$model$modules$m$weight
fitted$model$modules$m$bias

preds <- predict(fitted, valid_dl)
testPredictions = tibble(y = as_array(valid_ds$y), yhat = as_array(preds), x = as_array(valid_ds$x))
testPredictions

ggplot() +
  geom_point(aes(x = x, y = y, color = "y"), testPredictions, alpha = 0.1, shape = 16) +
  geom_point(aes(x = x, y = yhat, color = "yhat"), testPredictions, alpha = 0.1, shape = 16) +
  guides(color = "none")

#m
#x = torch_reshape(torch_tensor(1:4, dtype = torch_float32()), c(4, 1))
#x %>% m
#batch <- dataloader_make_iter(train_dl) %>% dataloader_next()
#model(batch$x)
#foo = torch_reshape(torch_tensor(seq(1:10), dtype = torch_double()), c(10, 1))
#bar = torch_reshape(torch_tensor(seq(2:11), dtype = torch_double()), c(10, 1))
#nnf_mse_loss(foo, bar)




# sanity checks
cuda_is_available()
cudaDeviceIndex = cuda_current_device()
torch_device("cuda", index = cudaDeviceIndex)

torch_tensor(na.omit(penguins[ , c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g", "year")]) %>% as.matrix())

torch_install_path()

# https://torch.mlverse.org/docs/
x <- torch_tensor(1, requires_grad = TRUE)
w <- torch_tensor(2, requires_grad = TRUE)
b <- torch_tensor(3, requires_grad = TRUE)
y <- w * x + b
y$backward()
x$grad
w$grad
b$grad

# guess the correlation https://torch.mlverse.org/start/guess_the_correlation/
torch_tensor(1)

train_indices <- 1:10000
val_indices <- 10001:15000
test_indices <- 15001:20000

add_channel_dim <- function(img) img$unsqueeze(1)
crop_axes <- function(img) transform_crop(img, top = 0, left = 21, height = 131, width = 130)

root <- file.path(tempdir(), "correlation")

train_ds <- guess_the_correlation_dataset(
  # where to unpack
  root = root,
  # additional preprocessing 
  transform = function(img) crop_axes(img) %>% add_channel_dim(),
  # don't take all data, but just the indices we pass in
  indexes = train_indices,
  download = TRUE
)

valid_ds <- guess_the_correlation_dataset(
  root = root,
  transform = function(img) crop_axes(img) %>% add_channel_dim(),
  indexes = val_indices,
  download = FALSE
)

test_ds <- guess_the_correlation_dataset(
  root = root,
  transform = function(img) crop_axes(img) %>% add_channel_dim(),
  indexes = test_indices,
  download = FALSE
)

length(train_ds)
length(valid_ds)
length(test_ds)

train_ds[1]

train_dl <- dataloader(train_ds, batch_size = 64, shuffle = TRUE)
length(train_dl)

batch <- dataloader_make_iter(train_dl) %>% dataloader_next()

dim(batch$x)
dim(batch$y)

par(mfrow = c(8,8), mar = rep(0, 4))

images <- as.array(batch$x$squeeze(2))

images %>%
  purrr::array_tree(1) %>%
  purrr::map(as.raster) %>%
  purrr::iwalk(~{plot(.x)})

batch$y %>% as.numeric() %>% round(digits = 2)

valid_dl <- dataloader(valid_ds, batch_size = 64)
length(valid_dl)

test_dl <- dataloader(test_ds, batch_size = 64)
length(test_dl)

torch_manual_seed(777)

net <- nn_module(
  
  "corr-cnn",
  
  initialize = function() {
    
    self$conv1 <- nn_conv2d(in_channels = 1, out_channels = 32, kernel_size = 3)
    self$conv2 <- nn_conv2d(in_channels = 32, out_channels = 64, kernel_size = 3)
    self$conv3 <- nn_conv2d(in_channels = 64, out_channels = 128, kernel_size = 3)
    
    self$fc1 <- nn_linear(in_features = 14 * 14 * 128, out_features = 128)
    self$fc2 <- nn_linear(in_features = 128, out_features = 1)
    
  },
  
  forward = function(x) {
    
    x %>% 
      self$conv1() %>%
      nnf_relu() %>%
      nnf_avg_pool2d(2) %>%
      
      self$conv2() %>%
      nnf_relu() %>%
      nnf_avg_pool2d(2) %>%
      
      self$conv3() %>%
      nnf_relu() %>%
      nnf_avg_pool2d(2) %>%
      
      torch_flatten(start_dim = 2) %>%
      self$fc1() %>%
      nnf_relu() %>%
      
      self$fc2()
  }
)

model <- net()
model(batch$x)


fitted <- net %>%
  setup(
    loss = function(y_hat, y_true) nnf_mse_loss(y_hat, y_true$unsqueeze(2)),
    optimizer = optim_adam
  )

# loss 0.0881 -> 0.0017
trainingStart = Sys.time()
fitted <- net %>%
  setup(
    loss = function(y_hat, y_true) nnf_mse_loss(y_hat, y_true$unsqueeze(2)),
    optimizer = optim_adam
  ) %>%
  fit(train_dl, epochs = 10, valid_data = test_dl)
Sys.time() - trainingStart

preds <- predict(fitted, test_dl)

preds <- preds$to(device = "cpu")$squeeze() %>% as.numeric()
test_dl <- dataloader(test_ds, batch_size = 5000)
targets <- (test_dl %>% dataloader_make_iter() %>% dataloader_next())$y %>% as.numeric()

df <- data.frame(preds = preds, targets = targets)

library(ggplot2)

ggplot(df, aes(x = targets, y = preds)) +
  geom_point(size = 0.1) +
  theme_classic() +
  xlab("true correlations") +
  ylab("model predictions")
